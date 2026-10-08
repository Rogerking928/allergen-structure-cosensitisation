"""晶片成分名 → UniProt accession（AlphaFold DB 以 accession 為鍵）。

權威來源是 WHO/IUIS Allergen Nomenclature（allergen.org）的逐一過敏原頁，
那裡每個同功異型都附 UniProt accession。UniProt 自己的全文搜尋不可靠：
多數 TrEMBL 條目的名稱欄位根本沒有 IUIS 名（Ara h 3 的 O82580 只叫 Glycinin），
反過來全文命中又常來自文獻標題。所以一律以 IUIS 為準，UniProt 只用來取序列。
"""
import json, re, time, requests
from io import StringIO
import pandas as pd
from common import load, ROOT

IUIS_LIST = "http://www.allergen.org/search.php?allergenname=&allergenname_exact=off&TaxSource=&TaxOrder=&foodallerg=all&bioname="
IUIS_VIEW = "http://www.allergen.org/viewallergen.php?aid={}"
CACHE = ROOT / "data" / "iuis"

# 多條鏈的過敏原：晶片上的成分是異源二聚體，兩條鏈都要留
MULTICHAIN = {"Fel_d_1": ["P30438", "P30440"]}
# 同一前驅物切出來的結構域，結構相似度會是人工的 1.0，之後要單獨處理
SUBCHAIN = {"Hev_b_6.01", "Hev_b_6.02", "Hev_b_6"}


def fetch(url, fname):
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / fname
    if not p.exists():
        for wait in [0, 5, 15, 45, 90]:          # allergen.org 會回 429，退避重試
            time.sleep(wait)
            r = requests.get(url, timeout=60, headers={"User-Agent": "allergen-structure-study/0.1"})
            if r.status_code == 200:
                break
        r.raise_for_status()
        p.write_text(r.text, encoding="utf-8")
        time.sleep(1.0)
    return p.read_text(encoding="utf-8")


def iuis_index():
    """IUIS 名 → aid。"""
    h = fetch(IUIS_LIST, "_index.html")
    return {n.strip(): a for a, n in re.findall(r'viewallergen\.php\?aid=(\d+)[^>]*>\s*([^<]+?)\s*<', h)}


def iuis_name(c):
    """Act_d_1 → 'Act d 1'；Cor_a_1.0101 → 'Cor a 1.0101'。不是單一成分則回 None。"""
    c = c.replace("_RUO", "")
    m = re.fullmatch(r"([A-Z][a-z]{2,3})_([a-z]{1,2})_(\d+(?:\.\d+)?)", c)
    return None if m is None else f"{m[1]} {m[2]} {m[3]}"


def isoallergens(aid):
    """該過敏原頁的同功異型表：[(iso_name, uniprot, pdb), ...]。"""
    h = fetch(IUIS_VIEW.format(aid), f"{aid}.html")
    try:
        tabs = pd.read_html(StringIO(h))
    except ValueError:
        return []
    out = []
    for t in tabs:
        if "Isoallergen and variants" in map(str, t.columns):
            for _, r in t.iterrows():
                # 這一欄的寫法很雜：'Q6PSU2-1'、'A5HII1 (variant V123F)'、空白分隔的兩個
                cell = str(r.get("UniProt", ""))
                accs = re.findall(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b", cell)
                out.append((str(r["Isoallergen and variants"]).strip(),
                            accs or [None],
                            str(r.get("PDB", "")).strip()))
    return out


def is_fragment(e):
    return "Fragment" in str(e.get("proteinDescription", {}).get("flag", ""))


def pick(isos, name):
    """選代表同功異型與 accession。

    晶片指名同功異型時優先用它；否則取 .0101。同一列可能列出多個 accession，
    其中常有只定序到一半的片段——片段會讓結構相似度整個失真，所以**完整序列優先**，
    再來才是 Swiss-Prot 與序列長度。
    """
    have = [i for i in isos if i[1] and i[1][0]]
    if not have:
        return None
    order = have
    if "." in name:
        order = [i for i in have if i[0] == name] or have
    else:
        order = [i for i in have if i[0].endswith(".0101")] or have
    cands = []
    for rank, (iso, accs, pdb) in enumerate(order + [h for h in have if h not in order]):
        for acc in accs:
            e = uniprot_entry(acc)
            if "sequence" not in e:        # 已刪除／合併的條目，查得到但沒有序列
                continue
            cands.append((rank, is_fragment(e), not e["entryType"].startswith("UniProtKB reviewed"),
                          -e["sequence"]["length"], iso, acc, pdb, e))
    if not cands:
        return None
    cands.sort(key=lambda t: t[:4])
    r = cands[0]
    return r[4], r[5], r[6], r[7]


UCACHE = ROOT / "data" / "uniprot"


def uniprot_entry(acc):
    """UniProt 條目，存成檔案快取——整輪要抓好幾百個，連線偶爾會被切斷，
    沒有快取的話每次重跑都得從頭來（2026-10-02 就這樣斷在最後幾個）。"""
    UCACHE.mkdir(parents=True, exist_ok=True)
    p = UCACHE / f"{acc}.json"
    if p.exists():
        return json.loads(p.read_text())
    last = None
    for wait in [0, 3, 10, 30]:
        time.sleep(wait)
        try:
            r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60)
            r.raise_for_status()
            p.write_text(r.text)
            time.sleep(0.1)
            return r.json()
        except Exception as e:
            last = e
    raise last


def main():
    _, _, comps = load()
    idx = iuis_index()
    rows = []
    for c in comps:
        name = iuis_name(c)
        if name is None:
            rows.append(dict(component=c, iuis=None, status="not_single_component"))
            continue
        accs, iso, pdb = [], None, ""
        if c in MULTICHAIN:
            accs, iso = MULTICHAIN[c], name
        else:
            aid = idx.get(name) or idx.get(name.split(".")[0])
            if aid is None:
                rows.append(dict(component=c, iuis=name, status="not_in_iuis"))
                continue
            got = pick(isoallergens(aid), name)
            if got is None:
                rows.append(dict(component=c, iuis=name, status="no_accession"))
                continue
            iso, accs, pdb = got[0], [got[1]], got[2]
        for acc in accs:
            e = uniprot_entry(acc)
            rows.append(dict(component=c, iuis=name, isoallergen=iso, status="ok", accession=acc,
                             reviewed=e["entryType"].startswith("UniProtKB reviewed"),
                             organism=e["organism"]["scientificName"],
                             length=e["sequence"]["length"], sequence=e["sequence"]["value"],
                             fragment="Fragment" in str(e.get("proteinDescription", {}).get("flag", "")),
                             pdb=pdb,
                             note="multichain" if c in MULTICHAIN else ("subchain" if c in SUBCHAIN else "")))

    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "out" / "component_uniprot.csv", index=False)
    u = df.drop_duplicates("component")
    print(u.status.value_counts().to_string())
    for s in ["not_in_iuis", "no_accession"]:
        bad = u[u.status == s]
        if len(bad):
            print(s, ":", bad.iuis.tolist())
    ok = df[df.status == "ok"]
    print("片段序列:", ok[ok.fragment].iuis.tolist())


if __name__ == "__main__":
    main()
