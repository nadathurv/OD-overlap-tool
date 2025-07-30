import pandas as pd
import logging
from pathlib import Path
from ..config import CDSCO_RAW, FDA_RAW

# Setup logging for data loading
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def load_raw():
    """Load raw CDSCO & FDA files (xlsx + csv)."""
    logging.info(f"Loading raw datasets: {CDSCO_RAW}, {FDA_RAW}")
    # Check raw files exist
    if not Path(CDSCO_RAW).exists():
        logging.error(f"CDSCO raw file not found: {CDSCO_RAW}")
        raise FileNotFoundError(f"CDSCO raw file not found: {CDSCO_RAW}")
    if not Path(FDA_RAW).exists():
        logging.error(f"FDA raw file not found: {FDA_RAW}")
        raise FileNotFoundError(f"FDA raw file not found: {FDA_RAW}")
    cdsco = pd.read_csv(CDSCO_RAW)
    fda = pd.read_csv(FDA_RAW)
    return cdsco, fda
