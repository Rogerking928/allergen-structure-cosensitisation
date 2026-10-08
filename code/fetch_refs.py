#!/usr/bin/env python3
"""Fetch and cache PubMed records for every reference cited in the manuscript.

Nothing in refs.py is written from memory: this script pulls the ESummary
record for each PMID and stores it in out/refs/verified.json, and refs.py is
generated from that cache by build_refs.py.

  python3 code/fetch_refs.py [extra_pmid ...]
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "refs")
CACHE = os.path.join(OUT, "verified.json")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def esummary(pmids):
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids), "retmode": "json"})
    with urllib.request.urlopen(url, timeout=60) as fh:
        return json.load(fh)["result"]


def main():
    # every PMID the manuscript cites is listed in build_refs.KEYS; nothing is
    # fetched from a study table here because this paper cites no included studies
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import build_refs
    pmids = list(build_refs.KEYS.values()) + sys.argv[1:]
    pmids = sorted(set(pmids))

    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo = [p for p in pmids if p not in cache]
    for i in range(0, len(todo), 50):
        chunk = todo[i:i + 50]
        res = esummary(chunk)
        for p in chunk:
            if p in res:
                cache[p] = res[p]
        time.sleep(0.4)
    json.dump(cache, open(CACHE, "w"), ensure_ascii=False, indent=1)
    print(f"cached {len(cache)} records in {CACHE}")
    for p in pmids:
        r = cache.get(p)
        print(f"{p}  {r['source'] if r else 'MISSING'}  {r['title'][:80] if r else ''}")


if __name__ == "__main__":
    main()
