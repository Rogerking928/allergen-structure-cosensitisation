"""抓 AlphaFold DB 預測結構（以 UniProt accession 為鍵），存 data/afdb/<acc>.pdb。

AFDB 的檔名帶版本號，所以先問 API 拿 URL，不要自己拼。
沒有模型的 accession 記在 out/afdb_missing.csv，之後要在流程圖裡交代。
"""
import time, requests
import pandas as pd
from common import ROOT

API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
OUT = ROOT / "data" / "afdb"


def fetch(acc):
    p = OUT / f"{acc}.pdb"
    if p.exists() and p.stat().st_size > 0:
        return "cached"
    r = requests.get(API.format(acc), timeout=60)
    if r.status_code != 200 or not r.json():
        return "no_model"
    meta = r.json()[0]
    url = meta.get("pdbUrl") or meta.get("cifUrl")
    s = requests.get(url, timeout=120)
    s.raise_for_status()
    p.write_bytes(s.content)
    time.sleep(0.2)
    return "downloaded"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(ROOT / "out" / "component_uniprot.csv")
    accs = sorted(m[m.status == "ok"].accession.dropna().unique())
    res = {}
    for a in accs:
        try:
            res[a] = fetch(a)
        except Exception as e:
            res[a] = f"error:{type(e).__name__}"
    s = pd.Series(res)
    print(s.value_counts().to_string())
    miss = s[s != "downloaded"][s != "cached"]
    miss.rename("reason").to_csv(ROOT / "out" / "afdb_missing.csv", index_label="accession")
    if len(miss):
        print("缺模型:", miss.index.tolist())


if __name__ == "__main__":
    main()
