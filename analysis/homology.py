"""成分對的同源關係，改用 UniProt 的 Pfam／SUPFAM 註解，不再靠手寫的家族關鍵字。

手寫家族標籤切得太細：7S 與 11S 都是 cupin、β-lactoglobulin 本身就是 lipocalin、
Asp f 3 與 Mala s 5 同為 peroxiredoxin——這些在舊版被算成「跨家族」，
造出了一批看似「結構像、家族不同」的假訊號（2026-10-02 逐對檢查發現）。

定義：兩個成分只要共享任一 Pfam family 或任一 SCOP superfamily（SUPFAM）即視為同源。
SUPFAM 對應的是超家族層級，能把 cupin、lipocalin 這類把多個 Pfam 收在一起的群併起來。
多鏈成分（Fel d 1）取兩條鏈註解的聯集。
"""
import json
import pandas as pd
from common import ROOT

CACHE = ROOT / "data" / "uniprot"
DBS = ("Pfam", "SUPFAM")


def annotations(acc):
    p = CACHE / f"{acc}.json"
    if not p.exists():
        return set()
    d = json.loads(p.read_text())
    return {f"{r['database']}:{r['id']}" for r in d.get("uniProtKBCrossReferences", [])
            if r["database"] in DBS}


def names(acc):
    p = CACHE / f"{acc}.json"
    if not p.exists():
        return {}
    d = json.loads(p.read_text())
    out = {}
    for r in d.get("uniProtKBCrossReferences", []):
        if r["database"] in DBS:
            nm = next((x["value"] for x in r.get("properties", []) if x["key"] == "EntryName"), "")
            out[f"{r['database']}:{r['id']}"] = nm
    return out


def table():
    m = pd.read_csv(ROOT / "out" / "component_uniprot.csv")
    m = m[(m.status == "ok") & m.accession.notna()]
    ann = {}
    for comp, g in m.groupby("component"):
        s = set()
        for acc in g.accession:
            s |= annotations(acc)
        ann[comp] = s
    return ann


def homologous(ann, a, b):
    A, B = ann.get(a, set()), ann.get(b, set())
    return bool(A & B)


if __name__ == "__main__":
    ann = table()
    n0 = sum(1 for v in ann.values() if not v)
    print(f"成分 {len(ann)}，沒有任何 Pfam/SUPFAM 註解 {n0}")
    for k, v in sorted(ann.items()):
        if not v:
            print("  無註解:", k)
