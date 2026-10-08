"""過敏原的蛋白家族標籤，取自 WHO/IUIS 頁面的 Biochemical name 欄（`data/iuis/` 快取）。

家族標籤只用來分組顯示與「同／跨家族」的分層，不參與估計；
原始的 Biochemical name 字串一併留在輸出裡，讀者能自行核對歸類。
"""
import re
from io import StringIO
import pandas as pd
from common import ROOT

CACHE = ROOT / "data" / "iuis"

# 關鍵字 → 家族標籤。順序有意義（先專一後一般）。
RULES = [
    (r"\bPR-10\b|Bet v 1 family|Bet v 1-like", "PR-10 (Bet v 1-like)"),
    (r"profilin", "Profilin"),
    (r"polcalcin|EF-hand.*2 EF|2 EF-hand", "Polcalcin (2 EF-hand)"),
    (r"parvalbumin", "Parvalbumin"),
    (r"tropomyosin", "Tropomyosin"),
    (r"lipid transfer|\bLTP\b|nsLTP", "nsLTP"),
    (r"2S (sulfur-rich )?(seed storage )?albumin|conglutin|napin", "2S albumin (prolamin)"),
    (r"group 5/21|group 5/21 allergen", "Mite group 5/21"),
    (r"antigen 5\b", "Antigen 5 (CAP superfamily)"),
    (r"oleosin", "Oleosin"),
    (r"7S|vicilin|cupin.*7S", "7S globulin (vicilin)"),
    (r"11S|legumin|glycinin|bicupin", "11S globulin (legumin)"),
    (r"lipocalin", "Lipocalin"),
    (r"serum albumin", "Serum albumin"),
    (r"pectate lyase", "Pectate lyase"),
    (r"beta-expansin|grass group (1|i)\b|group 1.*grass", "Grass group 1 (beta-expansin)"),
    (r"expansin|grass group (2|ii)|group 2.*grass", "Grass group 2/3 (expansin)"),
    (r"cysteine protease|papain", "Cysteine protease"),
    (r"serine protease|trypsin-like", "Serine protease"),
    (r"\bNPC2\b|ML domain|group 2 mite", "Mite group 2 (NPC2)"),
    (r"chitinase", "Chitinase"),
    (r"thaumatin|\bTLP\b", "Thaumatin-like (PR-5)"),
    (r"defensin", "Defensin"),
    (r"alpha-amylase.*inhibitor|trypsin inhibitor", "Alpha-amylase/trypsin inhibitor"),
    (r"gliadin|glutenin|prolamin", "Gluten prolamin"),
    (r"arginine kinase", "Arginine kinase"),
    (r"myosin light chain", "Myosin light chain"),
    (r"sarcoplasmic calcium", "Sarcoplasmic Ca-binding protein"),
    (r"troponin", "Troponin C"),
    (r"hevein", "Hevein-like"),
    (r"ole e 1|ole e 1-like|common olive group 1", "Ole e 1-like"),
    (r"lysozyme", "Lysozyme"),
    (r"transferrin|lactoferrin", "Transferrin"),
    (r"caseins?\b", "Casein"),
    (r"lactoglobulin", "Beta-lactoglobulin"),
    (r"lactalbumin", "Alpha-lactalbumin"),
    (r"ovomucoid|ovalbumin|ovotransferrin", "Egg white protein"),
    (r"uteroglobin|secretoglobin", "Secretoglobin"),
    (r"enolase", "Enolase"),
    (r"peroxisomal|thioredoxin|superoxide dismutase|cyclophilin", "Conserved cytosolic enzyme"),
]


def biochemical_names():
    """aid 頁 → {allergen name: biochemical name}。"""
    out = {}
    for p in CACHE.glob("*.html"):
        if p.stem.startswith("_"):
            continue
        h = p.read_text(encoding="utf-8", errors="replace")
        try:
            tabs = pd.read_html(StringIO(h))
        except ValueError:
            continue
        for t in tabs:
            if t.shape[1] == 2:
                d = dict(zip(t[0].astype(str).str.strip(), t[1].astype(str).str.strip()))
                name = d.get("Allergen name:")
                if name:
                    out[name] = d.get("Biochemical name:", "")
    return out


def family(bioname):
    s = (bioname or "").lower()
    for pat, lab in RULES:
        if re.search(pat, s, flags=re.I):
            return lab
    return "Other / unclassified"


def table():
    """component_uniprot.csv ＋ 家族標籤。"""
    m = pd.read_csv(ROOT / "out" / "component_uniprot.csv")
    bn = biochemical_names()
    m["biochemical_name"] = m.iuis.map(lambda x: bn.get(x, bn.get(str(x).split(".")[0], "")))
    m["family"] = m.biochemical_name.map(family)
    return m


if __name__ == "__main__":
    t = table()
    t.drop(columns=["sequence"], errors="ignore").to_csv(ROOT / "out" / "component_family.csv", index=False)
    ok = t[t.status == "ok"].drop_duplicates("component")
    print(ok.family.value_counts().to_string())
    print("\n未分類：")
    print(ok[ok.family == "Other / unclassified"][["iuis", "biochemical_name"]].to_string(index=False))
