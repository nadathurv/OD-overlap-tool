# Reduce to essential text normalization utilities
import re, unicodedata


def normalize(txt: str) -> str:
    if txt is None:
        return ""
    txt = unicodedata.normalize("NFKD", txt.lower())
    txt = re.sub(r"[^a-z0-9 ]", " ", txt)
    return re.sub(r"\s+", " ", txt).strip()


def jaccard(a: str, b: str) -> float:
    sa, sb = map(set, map(str.split, (normalize(a), normalize(b))))
    return len(sa & sb) / len(sa | sb | {" "})
