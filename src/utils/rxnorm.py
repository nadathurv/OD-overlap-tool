"""
Utilities for resolving drug names to RxNorm IDs using NLM RxNav API, with local caching.
"""
import json
import logging
import time
from pathlib import Path
import requests
from ..config import PROC

# Path for caching name->RxCUI mappings
tRXCUI_CACHE = PROC / "rxnorm_cache.json"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def _load_cache():
    if not tRXCUI_CACHE.exists():
        return {}
    try:
        return json.loads(tRXCUI_CACHE.read_text(encoding='utf-8'))
    except Exception as e:
        logging.error(f"Failed to load RxNorm cache: {e}")
        return {}


def _save_cache(cache: dict):
    try:
        tRXCUI_CACHE.parent.mkdir(exist_ok=True, parents=True)
        tRXCUI_CACHE.write_text(json.dumps(cache, indent=2), encoding='utf-8')
    except Exception as e:
        logging.error(f"Failed to save RxNorm cache: {e}")


def name_to_rxcui(name: str) -> str | None:
    """
    Resolve a drug name to a single RxCUI ID via the NLM RxNav API, with caching.

    Args:
        name: Drug name string (normalized)

    Returns:
        str or None: RxCUI identifier if found, else None
    """
    cache = _load_cache()
    if name in cache:
        return cache[name]

    rxcui = None
    try:
        url = f"https://rxnav.nlm.nih.gov/REST/rxcui.json?name={requests.utils.quote(name)}"
        resp = requests.get(url, timeout=5)
        data = resp.json()
        ids = data.get('idGroup', {}).get('rxnormId', [])
        if ids:
            rxcui = ids[0]
            logging.info(f"Mapped '{name}' to RxCUI {rxcui}")
    except Exception as e:
        logging.error(f"RxNorm lookup error for '{name}': {e}")
    cache[name] = rxcui
    _save_cache(cache)
    # Rate limiting
    time.sleep(0.1)
    return rxcui 