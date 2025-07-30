"""
Optimized fuzzy matching algorithms for drug name comparison.
"""

import pandas as pd
import textdistance
from rapidfuzz import fuzz
from typing import Dict, List, Tuple


# Optimized thresholds based on empirical analysis
OPTIMAL_JARO_THRESHOLD = 0.72      # Tightened for better precision
OPTIMAL_JACCARD_THRESHOLD = 0.1    # Optimal pre-filtering balance
OPTIMAL_TOKEN_THRESHOLD = 82       # RapidFuzz token-set ratio
OPTIMAL_RATIO_THRESHOLD = 82       # Levenshtein ratio


def jaro_winkler_similarity(s1: str, s2: str) -> float:
    """Calculate Jaro-Winkler similarity between two strings."""
    if not isinstance(s1, str) or not isinstance(s2, str):
        return 0.0
    return textdistance.jaro_winkler(s1, s2)


def is_high_confidence_match(jw: float, token: int, ratio: int,
                            jw_thresh: float = OPTIMAL_JARO_THRESHOLD,
                            token_thresh: int = OPTIMAL_TOKEN_THRESHOLD,
                            ratio_thresh: int = OPTIMAL_RATIO_THRESHOLD) -> bool:
    """
    Determine if a candidate pair is a strong match by requiring at least 
    two of three similarity metrics to meet or exceed their thresholds.
    
    Args:
        jw: Jaro-Winkler similarity score
        token: RapidFuzz token-set ratio score  
        ratio: Levenshtein ratio score
        jw_thresh: Jaro-Winkler threshold
        token_thresh: Token-set ratio threshold
        ratio_thresh: Ratio threshold
        
    Returns:
        bool: True if at least 2 of 3 metrics exceed thresholds
    """
    criteria = [
        jw >= jw_thresh,
        token >= token_thresh,
        ratio >= ratio_thresh,
    ]
    return sum(criteria) >= 2


def vectorized_jaccard_filter(query_name: str, candidate_names: pd.Series,
                            threshold: float = OPTIMAL_JACCARD_THRESHOLD) -> pd.Series:
    """
    Optimized Jaccard similarity filtering using vectorized operations.
    
    Args:
        query_name: Normalized query drug name
        candidate_names: Series of normalized candidate drug names
        threshold: Jaccard similarity threshold
        
    Returns:
        Boolean mask for candidates above threshold
    """
    query_tokens = set(query_name.split())
    
    def jaccard_score(candidate: str) -> float:
        if not candidate:
            return 0.0
        candidate_tokens = set(candidate.split())
        intersection = len(query_tokens & candidate_tokens)
        union = len(query_tokens | candidate_tokens)
        return intersection / union if union > 0 else 0.0
    
    # Vectorized application
    scores = candidate_names.apply(jaccard_score)
    return scores >= threshold


def create_match_record(cdsco_row: pd.Series, fda_row: pd.Series,
                       jw_score: float, token_score: int, ratio_score: int,
                       match_type: str = "Fuzzy") -> Dict:
    """
    Create a standardized match record dictionary.
    
    Args:
        cdsco_row: CDSCO dataset row
        fda_row: FDA dataset row  
        jw_score: Jaro-Winkler similarity score
        token_score: Token-set ratio score
        ratio_score: Levenshtein ratio score
        match_type: Type of match (RxNorm, Fuzzy, etc.)
        
    Returns:
        Dict containing standardized match data
    """
    return {
        "CDSCO Drug Name": cdsco_row["Drug Name"],
        "FDA Drug Name": fda_row["Drug Name"],
        "Similarity Score": jw_score,
        "Token Score": token_score,
        "Ratio Score": ratio_score,
        "Match Type": match_type,
        "CDSCO Approval Date": cdsco_row.get("Date of Approval", ""),
        "FDA Approval Date": fda_row.get("Date of Approval", ""),
        "CDSCO Indication": cdsco_row.get("Indication", ""),
        "FDA Indication": fda_row.get("Indication", ""),
    }


def optimized_fuzzy_matching(cdsco_df: pd.DataFrame, fda_df: pd.DataFrame,
                           synonyms: Dict[str, str] = None,
                           jw_threshold: float = OPTIMAL_JARO_THRESHOLD,
                           jaccard_threshold: float = OPTIMAL_JACCARD_THRESHOLD,
                           token_threshold: int = OPTIMAL_TOKEN_THRESHOLD,
                           ratio_threshold: int = OPTIMAL_RATIO_THRESHOLD) -> List[Dict]:
    """
    Optimized fuzzy matching algorithm with vectorized operations and early termination.
    
    Args:
        cdsco_df: CDSCO dataset with normalized names
        fda_df: FDA dataset with normalized names
        synonyms: Optional synonyms mapping
        jw_threshold: Jaro-Winkler threshold
        jaccard_threshold: Jaccard pre-filtering threshold
        token_threshold: Token-set ratio threshold
        ratio_threshold: Ratio threshold
        
    Returns:
        List of match dictionaries
    """
    if synonyms is None:
        synonyms = {}
    
    matches = []
    
    # Pre-compute FDA normalized names for vectorized operations
    fda_names = fda_df["drug_norm"].values
    
    for _, cdsco_row in cdsco_df.iterrows():
        # Apply synonym mapping
        query_name = cdsco_row["drug_norm"]
        query_name = synonyms.get(query_name, query_name)
        
        # Vectorized Jaccard pre-filtering
        jaccard_mask = vectorized_jaccard_filter(query_name, fda_df["drug_norm"], 
                                               jaccard_threshold)
        candidates = fda_df[jaccard_mask]
        
        # Early termination if no candidates pass Jaccard filter
        if candidates.empty:
            continue
            
        # Compute similarities for filtered candidates
        for _, fda_row in candidates.iterrows():
            jw_score = jaro_winkler_similarity(query_name, fda_row["drug_norm"])
            token_score = fuzz.token_set_ratio(query_name, fda_row["drug_norm"])
            ratio_score = fuzz.ratio(query_name, fda_row["drug_norm"])
            
            # Apply consensus matching
            if is_high_confidence_match(jw_score, token_score, ratio_score,
                                      jw_threshold, token_threshold, ratio_threshold):
                match_record = create_match_record(cdsco_row, fda_row, 
                                                 jw_score, token_score, ratio_score)
                matches.append(match_record)
    
    return matches


def process_id_matches(cdsco_df: pd.DataFrame, fda_df: pd.DataFrame,
                      id_column: str = "RxCUI") -> Tuple[List[Dict], pd.DataFrame]:
    """
    Process RxNorm ID-based matches and return remaining CDSCO entries.
    
    Args:
        cdsco_df: CDSCO dataset
        fda_df: FDA dataset
        id_column: Column name for ID matching (default: RxCUI)
        
    Returns:
        Tuple of (id_matches_list, remaining_cdsco_df)
    """
    id_matches = []
    
    if id_column in cdsco_df.columns and id_column in fda_df.columns:
        id_c = cdsco_df[cdsco_df[id_column].notna()]
        id_f = fda_df[fda_df[id_column].notna()]
        id_pairs = id_c.merge(id_f, on=id_column, suffixes=('_cdsco', '_fda'))
        
        for _, pair in id_pairs.iterrows():
            match_record = {
                'CDSCO Drug Name': pair['Drug Name_cdsco'],
                'FDA Drug Name': pair['Drug Name_fda'],
                'Similarity Score': 1.0,
                'Token Score': 100,
                'Ratio Score': 100,
                'Match Type': 'RxNorm',
                'CDSCO Approval Date': pair.get('Date of Approval_cdsco', ''),
                'FDA Approval Date': pair.get('Date of Approval_fda', ''),
                'CDSCO Indication': pair.get('Indication_cdsco', ''),
                'FDA Indication': pair.get('Indication_fda', ''),
            }
            id_matches.append(match_record)
        
        # Remove ID-matched entries from fuzzy matching pool
        matched_ids = set(id_pairs[id_column])
        remaining_cdsco = cdsco_df[~cdsco_df[id_column].isin(matched_ids)]
        
        return id_matches, remaining_cdsco
    
    return id_matches, cdsco_df 