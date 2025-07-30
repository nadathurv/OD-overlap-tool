# Orphan Drug Approval Overlap (CDSCO vs FDA)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

A high-performance pipeline that identifies overlapping drug approvals between India's CDSCO and the U.S. FDA orphan drug databases.

## How It Works: Technical Overview

### 1. **Data Cleaning & Standardization**
- **Input**: Raw CSV files from CDSCO and FDA
- **Process**: 
  - Normalize drug names (lowercase, remove special characters)
  - Strip dosage forms and strengths  
  - Standardize date formats
  - Handle combination drugs by splitting into individual APIs
- **Output**: Clean Parquet files for efficient processing

### 2. **Two-Phase Matching Algorithm**

#### Phase 1: RxNorm ID Matching (100% Accuracy)
- Direct match on RxNorm identifiers when available
- Guarantees perfect matches for drugs with standardized IDs

#### Phase 2: Fuzzy Name Matching
- **Pre-filtering**: Jaccard similarity (threshold: 0.1) for fast candidate selection
- **Consensus matching**: Requires 2 of 3 metrics to exceed thresholds:
  - Jaro-Winkler similarity ≥ 0.72
  - Token-set ratio ≥ 82
  - Levenshtein ratio ≥ 82
- **Synonym harmonization**: Maps drug variants to canonical names
- **1:1 correspondence**: Each drug matches at most once

### 3. **Performance Optimizations**
- Vectorized Jaccard filtering reduces O(n²) to O(n·m')
- Early termination when no candidates pass pre-filter
- Parquet format for 10x faster I/O than CSV

## Quick Start

### Prerequisites
- Python 3.8+
- pip

### Installation
```bash
git clone https://github.com/vnadathur/orphan-drug-overlap.git
cd orphan-drug-overlap
pip install -r requirements.txt
```

### Data Setup
Place raw CSV files in `data/raw/`:
- `cdsco_drugs.csv` - CDSCO drug list
- `fda_orphan_drugs.csv` - FDA orphan drug list

### Running the Pipeline

#### Standard Comparison
```bash
bash run_all.sh
```

#### With Combination Drug Splitting
```bash
bash run_all.sh --explode-combinations
```

#### Multi-Variant Analysis
```bash
bash run_multi_variant.sh
```

#### Custom Thresholds
```bash
python compare.py --jaro 0.75 --jaccard 0.15
```

## Output

The pipeline generates `data/processed/overlap.csv` with:
- **CDSCO Drug Name** & **FDA Drug Name**: Matched drug names
- **Similarity Score**: Jaro-Winkler score (0-1)
- **Token Score** & **Ratio Score**: Additional similarity metrics
- **Match Type**: "RxNorm" or "Fuzzy"
- **Approval Dates** & **Indications**: From both databases

## Advanced Usage

### Utility Scripts (in `utils/` directory)
- `generate_synonyms.py` - Create synonym mappings
- `check_synonyms.py` - Validate synonym coverage
- `sweep_thresholds.py` - Optimize threshold values

### API Reference
```python
from compare import run_comparison

# Run with custom parameters
results_df, output_path = run_comparison(
    cdsco_source="cdsco_clean",
    fda_source="fda_clean", 
    jaro=0.75,
    jaccard=0.15
)
```

## Performance Metrics
- Processes ~3,000 drugs in <5 seconds
- Achieves 95%+ precision with optimized thresholds
- Memory efficient with streaming Parquet operations

## License
MIT License - see [LICENSE](LICENSE) file