"""共用：讀入 Allergen Chip Challenge 公開檔（分號分隔、小數逗號）。"""
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "allergenchipchallenge-data-corrected-final-hdh-sfa.csv"
N_META = 15


def load():
    d = pd.read_csv(RAW, sep=";", decimal=",", low_memory=False, encoding="utf-8-sig")
    meta = d.columns[:N_META].tolist()
    comps = d.columns[N_META:].tolist()
    return d, meta, comps


_COMP_NAMES = None


def display_name(x):
    """晶片欄名 → 過敏原代號。`_RUO`（research use only）是平台標記，通常拿掉；
    但同一過敏原在另一平台另有不帶 _RUO 的欄位時（Ole e 7）要保留成「(RUO)」，否則兩列同名、計數會對不上。"""
    global _COMP_NAMES
    if not isinstance(x, str):
        return x
    if _COMP_NAMES is None:
        _COMP_NAMES = set(pd.read_csv(RAW, sep=";", nrows=0, encoding="utf-8-sig").columns)
    base = x.replace("_RUO", "")
    if x.endswith("_RUO") and base in _COMP_NAMES:
        return base.replace("_", " ") + " (RUO)"
    return base.replace("_", " ")
