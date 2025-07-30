from pathlib import Path

# Root directories
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
REP = ROOT / "reports"

# File names – adjust to your actual csv/xlsx names
CDSCO_RAW = RAW / "cdsco_drugs.csv"
FDA_RAW = RAW / "fda_orphan_drugs.csv"

CDSCO_CLEAN = PROC / "cdsco_clean.parquet"
FDA_CLEAN = PROC / "fda_clean.parquet"

# Added for combination drug support
CDSCO_EXPLODED = PROC / "cdsco_exploded.parquet"
FDA_VOCAB_PATH = PROC / "fda_api_vocab.pkl"
