"""
Utilities for splitting combination drug names into individual API components.
"""

import re
import logging
import textdistance
from typing import List, Set, Tuple
from .text import normalize

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

# Constants for matching thresholds
VOCAB_MATCH_THRESHOLD = 0.9  # Threshold for considering a split valid when matching against vocabulary
FINAL_FILTER_THRESHOLD = 0.8  # Threshold for filtering out fragments that don't match any vocabulary term

# Special case mappings for known combinations that have alternative canonical names
SPECIAL_COMBINATIONS = {
    "folinic acid and calcium": "calcium folinate",
    "folinic acid with calcium": "calcium folinate",
    "folinic acid & calcium": "calcium folinate",
    "folinic acid + calcium": "calcium folinate",
    "folinic acid, calcium": "calcium folinate",
    "leucovorin": "calcium folinate",  # Leucovorin is another name for calcium folinate
    "calcium folinate": "calcium folinate",  # Include the canonical form itself
}


def _extract_parentheticals(text: str) -> Tuple[str, List[str]]:
    """
    Extract top-level parenthetical content from text using a stack-based parser.
    
    Args:
        text: Raw drug name possibly containing parenthetical content
        
    Returns:
        Tuple containing:
            - Modified text with parentheticals removed
            - List of extracted parenthetical content
    """
    if not text:
        return "", []
    
    result = []
    stack = []
    start_idx = -1
    backbone = list(text)
    
    for i, char in enumerate(text):
        if char == '(':
            if not stack:  # Start of a top-level parenthesis
                start_idx = i
            stack.append(i)
        elif char == ')' and stack:
            stack.pop()
            if not stack:  # End of a top-level parenthesis
                # Extract the content without the parentheses
                content = text[start_idx + 1:i].strip()
                if content:
                    result.append(content)
                # Mark for removal in the backbone
                for j in range(start_idx, i + 1):
                    backbone[j] = ' '
    
    # Rebuild the backbone string without parenthetical content
    backbone_text = ''.join(backbone).strip()
    backbone_text = re.sub(r'\s+', ' ', backbone_text)
    
    return backbone_text, result


def _split_by_delimiters(text: str) -> List[str]:
    """
    Split text by common drug delimiter characters (comma, slash, plus, ampersand).
    
    Args:
        text: Text to split
        
    Returns:
        List of split components
    """
    if not text:
        return []
    
    # Replace common delimiters with a standard delimiter
    standardized = text
    standardized = re.sub(r'\s*[,/+&]\s*', '|', standardized)
    standardized = re.sub(r'\s+and\s+', '|', standardized)
    standardized = re.sub(r'\s+with\s+', '|', standardized)  # Add "with" as a delimiter
    
    # Split on the standard delimiter
    parts = [part.strip() for part in standardized.split('|')]
    return [part for part in parts if part]


def _best_vocab_match(term: str, vocab: Set[str]) -> Tuple[str, float]:
    """
    Find the best matching vocabulary term for a given string.
    
    Args:
        term: Term to match
        vocab: Set of vocabulary terms to match against
        
    Returns:
        Tuple containing:
            - Best matching vocabulary term
            - Jaro-Winkler similarity score
    """
    if not term or not vocab:
        return "", 0.0
    
    norm_term = normalize(term)
    
    # Check for special combinations first
    if norm_term in SPECIAL_COMBINATIONS and SPECIAL_COMBINATIONS[norm_term] in vocab:
        return SPECIAL_COMBINATIONS[norm_term], 1.0
    
    best_score = 0.0
    best_match = ""
    
    for vocab_term in vocab:
        score = textdistance.jaro_winkler(norm_term, vocab_term)
        if score > best_score:
            best_score = score
            best_match = vocab_term
            
    return best_match, best_score


def _handle_and_splits(parts: List[str], vocab: Set[str]) -> List[str]:
    """
    Handle special case of 'and' within terms by only splitting when both sides match vocabulary.
    
    Args:
        parts: List of potential drug components, some of which might contain 'and'
        vocab: Set of vocabulary terms to match against
        
    Returns:
        List of drug components with appropriate 'and' splitting
    """
    result = []
    
    for part in parts:
        # First check for special combinations
        norm_part = normalize(part)
        if norm_part in SPECIAL_COMBINATIONS and SPECIAL_COMBINATIONS[norm_part] in vocab:
            result.append(part)  # Keep intact, will be mapped in final normalization
            continue
            
        # Check if this part contains ' and ' that wasn't already split
        if ' and ' in part.lower():
            # First check if the whole string closely matches a vocabulary term
            whole_term_match, whole_term_score = _best_vocab_match(part, vocab)
            if whole_term_score >= VOCAB_MATCH_THRESHOLD:
                # The whole term is a known API, keep it intact
                result.append(part)
                continue
                
            # Try splitting on 'and'
            potential_split = re.split(r'\s+and\s+', part, flags=re.IGNORECASE)
            
            if len(potential_split) == 2:  # Simple "X and Y" case
                left, right = potential_split
                left_match, left_score = _best_vocab_match(left, vocab)
                right_match, right_score = _best_vocab_match(right, vocab)
                
                # Only split if both sides match vocabulary terms well
                if left_score >= VOCAB_MATCH_THRESHOLD and right_score >= VOCAB_MATCH_THRESHOLD:
                    result.extend([left, right])
                else:
                    result.append(part)  # Keep as is - "and" is likely part of a single API
            else:
                result.append(part)  # Complex case, keep as is
        else:
            result.append(part)
            
    return result


def split_apis(drug_name: str, vocab: Set[str]) -> List[str]:
    """
    Split a potentially combination drug name into individual API components.
    
    Args:
        drug_name: Raw drug name, possibly containing multiple APIs
        vocab: Set of normalized drug names to validate against
        
    Returns:
        List of individual API components
    """
    if not drug_name or not vocab:
        return []
    
    # Check if the entire drug name is a special combination first
    norm_drug = normalize(drug_name)
    if norm_drug in SPECIAL_COMBINATIONS and SPECIAL_COMBINATIONS[norm_drug] in vocab:
        canonical = SPECIAL_COMBINATIONS[norm_drug].title()
        return [canonical]
    
    # Quick exit for simple drug names without delimiters or parentheses
    if '(' not in drug_name and ',' not in drug_name and '/' not in drug_name and '+' not in drug_name and '&' not in drug_name and ' and ' not in drug_name.lower() and ' with ' not in drug_name.lower():
        # No indication of combination, return as single component
        return [drug_name.title()]
    
    # Special case for nested parentheses like "ABC/3TC (Abacavir/Lamivudine)"
    # Look for drug names directly in the string
    components = []
    for vocab_term in vocab:
        if len(vocab_term) > 4:  # Skip very short terms to avoid false matches
            # Create a pattern that matches the vocab term with word boundaries
            pattern = r'\b' + re.escape(vocab_term) + r'\b'
            if re.search(pattern, norm_drug, re.IGNORECASE):
                components.append(vocab_term.title())
    
    # If we found likely drug names embedded in the string, return them
    if len(components) >= 2:
        # Remove duplicates
        unique_components = []
        for comp in components:
            if comp not in unique_components:
                unique_components.append(comp)
        return unique_components
    
    # Pre-process complex combinations by looking for known special combinations
    # within the larger string
    for pattern, replacement in SPECIAL_COMBINATIONS.items():
        # Check if the pattern is in the drug name
        pattern_regex = r'\b' + re.escape(pattern) + r'\b'
        if re.search(pattern_regex, norm_drug, re.IGNORECASE):
            # Replace the special pattern with a placeholder
            norm_drug = re.sub(pattern_regex, "___PLACEHOLDER___", norm_drug, flags=re.IGNORECASE)
            drug_name = re.sub(pattern_regex, "___PLACEHOLDER___", drug_name, flags=re.IGNORECASE)
    
    # Step A: Extract parenthetical content
    backbone, parentheticals = _extract_parentheticals(drug_name)
    
    # Process both backbone and parentheticals
    components = []
    
    # Step B: Split backbone on unambiguous delimiters
    if backbone:
        backbone_parts = _split_by_delimiters(backbone)
        # Replace placeholders with the original special combination
        processed_parts = []
        for part in backbone_parts:
            if "___PLACEHOLDER___" in part:
                # Find the matching special combination
                for pattern, replacement in SPECIAL_COMBINATIONS.items():
                    if replacement in vocab:
                        processed_parts.append(pattern)
                        break
            else:
                processed_parts.append(part)
        components.extend(processed_parts)
    
    # Process parentheticals - only keep those that match vocabulary terms well
    # and ignore likely abbreviations (short strings like ABC, 3TC, etc.)
    for p in parentheticals:
        # Skip likely abbreviations (short all-caps strings or strings with numbers)
        if (len(p) <= 3) or (p.isupper() and len(p) <= 5) or any(c.isdigit() for c in p):
            continue
            
        _, score = _best_vocab_match(p, vocab)
        if score >= VOCAB_MATCH_THRESHOLD:
            # This parenthetical content looks like an API, keep it
            components.append(p)
        else:
            # This might be a compound in parentheses, try to split it
            p_parts = _split_by_delimiters(p)
            if p_parts:
                components.extend(p_parts)
            else:
                # Check if there are drug names directly in the parenthetical
                norm_p = normalize(p)
                p_components = []
                for vocab_term in vocab:
                    if len(vocab_term) > 4:
                        pattern = r'\b' + re.escape(vocab_term) + r'\b'
                        if re.search(pattern, norm_p, re.IGNORECASE):
                            p_components.append(vocab_term.title())
                if p_components:
                    components.extend(p_components)
    
    # Step C: Handle special case of 'and' within terms
    components = _handle_and_splits(components, vocab)
    
    # Step D: Normalize each component
    normalized_components = []
    for comp in components:
        # First normalize to lowercase for matching
        norm = normalize(comp)
        
        # Check for special combinations
        if norm in SPECIAL_COMBINATIONS and SPECIAL_COMBINATIONS[norm] in vocab:
            canonical = SPECIAL_COMBINATIONS[norm].title()
            normalized_components.append(canonical)
            continue
        
        # Then get the best title-cased version from vocabulary or use title-cased original
        if norm:
            best_match, score = _best_vocab_match(norm, vocab)
            # Step E: Filter out components that don't match vocabulary well
            if score >= FINAL_FILTER_THRESHOLD:
                # Use either the matched vocab term or a title-cased version of the original
                if score >= VOCAB_MATCH_THRESHOLD:
                    normalized_components.append(best_match.title())
                else:
                    normalized_components.append(comp.title())
    
    # Remove duplicates while preserving order
    unique_components = []
    for comp in normalized_components:
        if comp not in unique_components:
            unique_components.append(comp)
    
    return unique_components 