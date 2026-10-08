#!/usr/bin/env python3
"""Generate code/refs.py from the verified PubMed cache.

Entries are written by machine from out/refs/verified.json and the ESummary and
EFetch records behind it, so no citation in the manuscript is typed from
memory. Items that have no PubMed record (a guideline document, an article
indexed nowhere, the registries themselves) are held in MANUAL below, each with
the source that was checked.

  python3 code/fetch_refs.py && python3 code/build_refs.py
"""
import json
import os
import re
import urllib.parse
import urllib.request
from xml.etree import ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "out", "refs", "verified.json")
XMLCACHE = os.path.join(ROOT, "out", "refs", "efetch.xml")
DEST = os.path.join(ROOT, "code", "refs.py")

# citekey -> PMID. The key is author-year, which is what the manuscript source
# writes between braces; the number in the printed list comes from the order of
# first appearance and is assigned by number_refs.py.
KEYS = {
    # the data and the field it sits in
    "martinroche2026": "41016485",   # Allergen Chip Challenge open database, JACI 2026
    "dramburg2023": "37186333",      # EAACI Molecular Allergology User's Guide 2.0
    "radauer2008": "18395549",       # allergens fall into few protein families
    "aalberse2000": "10932064",      # structural biology of allergens
    # co-sensitisation already mapped on multiplex arrays
    "scala2010": "20214663",         # ISAC, 23,077 patients, clustering by homology and source
    "scala2011": "21949785",         # profilin / PR-10 / tropomyosin, 3,113 patients
    "palacin2012": "23272072",       # co-sensitisation graph of LTPs
    "zou2026": "42333264",           # Ising co-sensitisation network (extracts, children)
    # structure resources and methods
    "negi2023": "37781674",          # SDAP 2.0, AlphaFold models for allergens
    "varadi2024": "37933859",        # AlphaFold DB 2024
    "vankempen2024": "37156916",     # Foldseek
    "zhang2004": "15476259",         # TM-score
    "pomes2018": "29625844",         # WHO/IUIS allergen nomenclature
    "blum2025": "39565202",          # InterPro 2025
    "paysan2025": "39540428",        # Pfam 2025
    "jumper2021": "34265844",        # AlphaFold
    "kabasser2024": "36331131",      # co-sensitisation to non-homologous allergens of one source
    "zhang2005": "15849316",         # TM-align 演算法本身（本文全程用它）
    "andreeva2020": "31724711",      # SCOP 2020（同源判定的另一半來源）
    # 平台與閾值
    "heffler2018": "29743964",       # ALEX macroarray
    "lupinek2014": "24161540",       # ISAC / MeDALL allergen-chip
    # 那 116 對裡已被實驗描述過的關係（審閱要求逐一引用，避免讀成新發現）
    "hilger2012": "22791068",        # 動物 lipocalin 過敏原
    "pichler2015": "25978036",       # 果膠裂解酶花粉過敏原的交叉反應（Cry j 1 / Cup a 1 / Jun a 1）
    "walnut2024": "39083591",        # 核桃—榛果 2S albumin 交叉反應
    "cupin2005": "15935274",         # cupin 超家族的同三聚體過敏原（7S/11S）
}

MANUAL = {
    # The data release itself. Checked on data.gouv.fr 2026-10-03: producer Société
    # Française d'Allergologie, last updated 26 July 2024, Licence Ouverte 2.0.
    "datagouv2024":
        "Société Française d'Allergologie. (2024). Allergen Chip Challenge "
        "[dataset]. data.gouv.fr; last updated 26 July 2024; Licence Ouverte / "
        "Open Licence 2.0. https://www.data.gouv.fr/datasets/allergen-chip-challenge",
    # The sequence criterion. Codex Alimentarius guideline CAC/GL 45-2003, Annex 1
    # (assessment of possible allergenicity): >35% identity over a window of 80 amino
    # acids. Wording checked against secondary summaries (AllergenOnline help page).
    "codex2003":
        "Codex Alimentarius Commission. (2003). Guideline for the conduct of food "
        "safety assessment of foods derived from recombinant-DNA plants (CAC/GL "
        "45-2003), Annex 1: Assessment of possible allergenicity. Rome: FAO/WHO.",
}

# Full journal titles as the publishers print them; the NLM record carries a
# parenthetical or an "official publication of ..." tail for these.
JOURNAL = {
    "BMJ (Clinical research ed.)": "BMJ",
    "The journal of allergy and clinical immunology. In practice":
        "Journal of Allergy and Clinical Immunology: In Practice",
    "Annals of allergy, asthma & immunology : official publication of the "
    "American College of Allergy, Asthma, & Immunology":
        "Annals of Allergy, Asthma & Immunology",
    "Jornal brasileiro de pneumologia : publicacao oficial da Sociedade "
    "Brasileira de Pneumologia e Tisilogia": "Jornal Brasileiro de Pneumologia",
    "Allergology international : official journal of the Japanese Society of "
    "Allergology": "Allergology International",
    "Postepy dermatologii i alergologii":
        "Postepy Dermatologii i Alergologii",
    "Revista portuguesa de pneumologia": "Revista Portuguesa de Pneumologia",
}

SMALL = {"a", "an", "and", "as", "at", "but", "by", "for", "in", "of", "on",
         "or", "the", "to", "with", "from"}


def journal_title(raw):
    """Full journal title in headline case, as the house style prints it."""
    raw = " ".join(raw.split())
    if raw in JOURNAL:
        return JOURNAL[raw]
    raw = re.sub(r"^The\s+", "", raw).rstrip(".")
    words = raw.split()
    out = []
    for i, w in enumerate(words):
        if w.isupper() or any(c.isdigit() for c in w):
            out.append(w)
        elif i and w.lower() in SMALL:
            out.append(w.lower())
        else:
            out.append(w[0].upper() + w[1:])
    return " ".join(out)


def expand_pages(pages):
    """PubMed abbreviates the closing page: 974-8 is pages 974 to 978."""
    m = re.match(r"^(\d+)-(\d+)$", pages)
    if m and len(m.group(2)) < len(m.group(1)):
        start, end = m.groups()
        end = start[:len(start) - len(end)] + end
        return f"{start}\u2013{end}"
    return pages.replace("-", "\u2013")


def author(name):
    """'Plouin PF' -> 'Plouin, P. F.'"""
    parts = name.rsplit(" ", 1)
    if len(parts) == 1:
        return name
    surname, initials = parts
    if not initials.isupper() or len(initials) > 4:
        return name
    return surname + ", " + " ".join(c + "." for c in initials)


def fetch_xml(pmids):
    if os.path.exists(XMLCACHE):
        return ET.parse(XMLCACHE).getroot()
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids),
                                     "retmode": "xml"}))
    with urllib.request.urlopen(url, timeout=120) as fh:
        data = fh.read()
    open(XMLCACHE, "wb").write(data)
    return ET.fromstring(data)


def main():
    cache = json.load(open(CACHE))
    pmids = sorted(set(KEYS.values()))
    root = fetch_xml(pmids)
    full = {}
    for art in root.iter("PubmedArticle"):
        pmid = art.findtext(".//PMID")
        full[pmid] = art.findtext(".//Journal/Title") or ""

    entries = {}
    for key, pmid in KEYS.items():
        r = cache[pmid]
        names = [author(a["name"]) for a in r["authors"]
                 if a.get("authtype") == "Author"]
        if len(names) > 6:
            who = ", ".join(names[:6]) + ", et al."
        elif len(names) > 1:
            who = ", ".join(names[:-1]) + ", & " + names[-1]
        else:
            who = names[0] if names else ""
        year = re.match(r"(\d{4})", r["pubdate"]).group(1)
        title = r["title"].rstrip(".")
        jour = journal_title(full.get(pmid, r["source"]))
        vol = r.get("volume", "")
        iss = r.get("issue", "")
        pages = expand_pages(r.get("pages") or "")
        doi = ""
        for aid in r.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid["value"]
        where = jour
        if vol:
            where += f", {vol}"
            if iss:
                where += f"({iss})"
        if pages:
            where += f", {pages}"
        stop = "" if title.endswith(("?", "!")) else "."
        txt = f"{who} ({year}). {title}{stop} {where}."
        if doi:
            txt += f" https://doi.org/{doi}"
        entries[key] = txt

    entries.update(MANUAL)

    with open(DEST, "w") as fh:
        fh.write('"""Reference pool for the allergen structure and co-sensitisation brief report.\n\n'
                 "GENERATED by code/build_refs.py from out/refs/verified.json, which\n"
                 "holds the PubMed ESummary and EFetch records retrieved for every cited\n"
                 "PMID. Do not edit by hand: add the PMID to KEYS in build_refs.py and\n"
                 "regenerate, or add an item with no PubMed record to MANUAL there.\n\n"
                 "Numbers in the printed list are assigned by order of first appearance\n"
                 'in the manuscript, by number_refs.py."""\n\n')
        fh.write("POOL = {\n")
        for key in sorted(entries):
            body = entries[key].replace('"', '\\"')
            fh.write(f'    "{key}":\n')
            line = ""
            for word in body.split(" "):
                if len(line) + len(word) > 66:
                    fh.write(f'        "{line}"\n')
                    line = ""
                line += word + " "
            fh.write(f'        "{line.rstrip()}",\n')
        fh.write("}\n")
    print(f"wrote {DEST} with {len(entries)} entries")


if __name__ == "__main__":
    main()
