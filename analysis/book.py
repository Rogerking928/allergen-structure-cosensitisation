"""圖表本（Word）＋ 700 dpi PNG ＋ Source_data.xlsx，全部放同一個交付資料夾。

版式照嚴重氣喘／Bcc 那兩篇量到的規格：
Letter 橫向、左右 1.0／上下 0.9 吋（照 Bcc 那篇）、每項新起一頁；標題 TNR 11.5 pt 粗體、單行距、後距 4 pt；
表格滿版 8 pt、表頭粗體、三條黑線（single、sz 8、000000）、無直線；
表註 9 pt、圖說 9.5 pt，正體；圖：標題 → 9 pt 斜體圖說 → 置中圖片。先全部表、再全部圖。
**只交 PNG，不產 PDF。**

用法：python3 book.py [輸出資料夾]
"""
import os
import json, os, re, shutil, sys
import numpy as np
import pandas as pd
from PIL import Image
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Inches, Mm
from openpyxl import Workbook
from openpyxl.styles import Font as XFont, Alignment, PatternFill

from common import load, ROOT
import families

OUT = ROOT / "out"
FONT = "Times New Roman"
# 照 Bcc 那篇（life_disruption/docbuild.py）：Letter 橫向、左右 1.0 吋、上下 0.9 吋
PAGE_W, PAGE_H = Inches(11.0), Inches(8.5)
MARGIN_LR, MARGIN_TB = Inches(1.0), Inches(0.9)
BOX_W, BOX_H = 8.2, 5.4
MAX_SCALE = 1.5
# 交付資料夾：環境變數 DELIVERABLE_DIR 優先，其次命令列第一個參數，最後用本機預設
DEFAULT_DEST = os.environ.get(
    "DELIVERABLE_DIR", "/mnt/c/Users/roger/Desktop/還沒動/過敏原結構_共同致敏")

# 物種學名要斜體；過敏原代號（Bet v 1）照 IUIS 慣例用正體
ITALIC_RE = re.compile(r"\b([A-Z][a-z]+ [a-z]{3,}(?: [a-z]{3,})?)\b")


# ------------------------------------------------------------------ docx
def run(p, text, size=None, bold=None, italic=None):
    r = p.add_run(text)
    r.font.name = FONT
    r._element.get_or_add_rPr().get_or_add_rFonts().set(qn("w:eastAsia"), FONT)
    if size: r.font.size = Pt(size)
    if bold is not None: r.bold = bold
    if italic is not None: r.italic = italic
    return r


def run_sci(p, text, size, bold=None):
    """把學名切出來設斜體，其餘正體。"""
    for part in re.split(r"(<i>.*?</i>)", text):
        if part.startswith("<i>"):
            run(p, part[3:-4], size, bold=bold, italic=True)
        elif part:
            run(p, part, size, bold=bold, italic=False)


def para(doc, before=0, after=0, align=None, page_break=False, keep=False):
    p = doc.add_paragraph()
    f = p.paragraph_format
    f.space_before, f.space_after, f.line_spacing = Pt(before), Pt(after), 1.0
    f.page_break_before, f.keep_with_next, f.keep_together = page_break, keep, keep
    if align is not None: p.alignment = align
    return p


def borders(cell, *edges):
    tcPr = cell._tc.get_or_add_tcPr()
    for old in tcPr.findall(qn("w:tcBorders")): tcPr.remove(old)
    b = OxmlElement("w:tcBorders")
    for e in edges:
        el = OxmlElement(f"w:{e}")
        for k, v in (("val", "single"), ("sz", "8"), ("space", "0"), ("color", "000000")):
            el.set(qn(f"w:{k}"), v)
        b.append(el)
    tcPr.append(b)


def add_table(doc, item, first):
    p = para(doc, after=4, page_break=not first, keep=True)
    label, rest = item["caption"].split(". ", 1)
    run(p, label + ". ", 11.5, bold=True); run(p, rest, 11.5, bold=True)
    rows = item["rows"]
    ncol = len(rows[0])
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    first_w = min(3.6, BOX_W * 0.34) if ncol > 2 else BOX_W * 0.55
    widths = [first_w] + [(BOX_W - first_w) / (ncol - 1)] * (ncol - 1)
    for i, row in enumerate(rows):
        is_group = i > 0 and all(c == "" for c in row[1:])
        for j, text in enumerate(row):
            c = t.cell(i, j)
            c.width = Inches(widths[j])
            cp = c.paragraphs[0]
            cp.paragraph_format.space_after = Pt(0); cp.paragraph_format.line_spacing = 1.0
            cp.alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER
            run_sci(cp, str(text), 8, bold=(i == 0 or (is_group and j == 0)))
    last = len(rows) - 1
    for j in range(ncol):
        borders(t.cell(0, j), "top", "bottom")
        if last >= 1: borders(t.cell(1, j), *(("top", "bottom") if last == 1 else ("top",)))
        if last > 1: borders(t.cell(last, j), "bottom")
    trPr = t.rows[0]._tr.get_or_add_trPr()
    h = OxmlElement("w:tblHeader"); h.set(qn("w:val"), "true"); trPr.append(h)
    if item.get("note"):
        n = para(doc, before=3)
        run_sci(n, item["note"], 9)


def add_figure(doc, item, first):
    p = para(doc, after=4, page_break=not first, keep=True)
    run(p, item["caption"], 11.5, bold=True)
    if item.get("note"):
        n = para(doc, after=6, keep=True)
        run_sci(n, item["note"], 9.5)
    w_px, h_px = Image.open(item["path"]).size
    w_in, h_in = w_px / 700, h_px / 700
    s = min(MAX_SCALE, BOX_W / w_in, BOX_H / h_in)
    ip = para(doc, align=WD_ALIGN_PARAGRAPH.CENTER)
    ip.add_run().add_picture(str(item["path"]), width=Inches(w_in * s), height=Inches(h_in * s))


def build(items, path, supp=()):
    """一份 Word：主文的表、主文的圖，接著補充的表、補充的圖（使用者 2026-10-04：圖表 Word 一個就好）。"""
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name, st.font.size = FONT, Pt(8)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = PAGE_W, PAGE_H
    sec.left_margin = sec.right_margin = MARGIN_LR
    sec.top_margin = sec.bottom_margin = MARGIN_TB
    tabs = [i for i in items if i["kind"] == "table"]
    figs = [i for i in items if i["kind"] == "figure"]
    stabs = [i for i in supp if i["kind"] == "table"]
    sfigs = [i for i in supp if i["kind"] == "figure"]
    for k, it in enumerate(tabs + figs + stabs + sfigs):
        (add_table if it["kind"] == "table" else add_figure)(doc, it, k == 0)
    doc.save(path)
    return len(tabs) + len(stabs), len(figs) + len(sfigs)


# ------------------------------------------------------------------ 內容
def disp(name):
    """晶片欄名 → 過敏原代號。`_RUO`（research use only）是平台的標記，不是名字的一部分。"""
    from common import display_name
    return display_name(str(name))


def _or_ci(o, lo, hi):
    """勝算比與區間：低於 1 時一位小數會變成 0.1 (0.0 to 0.4)，改兩位。"""
    nd = 2 if hi < 1 else 1
    return f"{o:.{nd}f} ({lo:.{nd}f} to {hi:.{nd}f})"


def fmt(x, nd=2):
    return "-" if pd.isna(x) else f"{x:.{nd}f}"


def pairs():
    import figs
    return figs.load_pairs()


def _n(v):
    """計數的分位數：整數就寫整數，內插出 .5 就照寫（28.5 若寫成 28，和百分比那列對不上）。"""
    return f"{v:.0f}" if float(v).is_integer() else f"{v:.1f}"


def table1():
    d, meta, comps = load()
    rows = [["Characteristic", "ISAC v1", "ISAC v2", "ALEX2", "All"]]
    import figstyle as F
    cols = F.CHIP_ORDER
    def keeps(g):
        x = g[comps].replace(-1, np.nan)
        return [c for c in comps if x[c].notna().sum() >= 0.8 * len(g)]
    union = sorted({c for ch in cols for c in keeps(d[d.Chip_Type == ch])})

    def per_patient(g):
        x = g[comps].replace(-1, np.nan)
        keep = keeps(g)
        b = (x[keep] >= 0.3).sum(axis=1)
        return b, 100 * b / len(keep)

    def col(g, ncomp):
        # 「All」欄要把各平台自己的值合併（各平台用自己測得的成分當分母），
        # 不可拿全體重算 keeps——那只剩三平台共同的成分，曾算出 81.4%、8 (2-16) 的錯值
        parts = [per_patient(g[g.Chip_Type == ch]) for ch in cols if (g.Chip_Type == ch).any()]
        b = pd.concat([p_[0] for p_ in parts])
        pct = pd.concat([p_[1] for p_ in parts])
        age = g.Age.dropna()
        return dict(n=f"{len(g):,}", comps=f"{ncomp}",
                    age=f"{age.median():.0f} ({age.quantile(.25):.0f}-{age.quantile(.75):.0f})",
                    sens=f"{(g.Sensitization == 1).mean() * 100:.1f}",
                    gender1=f"{(g.Gender == 1).mean() * 100:.1f}",
                    pos1=f"{(b >= 1).mean() * 100:.1f}",
                    breadth=f"{_n(b.median())} ({_n(b.quantile(.25))}-{_n(b.quantile(.75))})",
                    pctpos=f"{pct.median():.1f} ({pct.quantile(.25):.1f}-{pct.quantile(.75):.1f})")
    cs = {c: col(d[d.Chip_Type == c], len(keeps(d[d.Chip_Type == c]))) for c in cols}
    cs["All"] = col(d, len(union))

    def r(label, key):
        rows.append([label] + [cs[c][key] for c in cols + ["All"]])
    r("Patients, n", "n")
    r("Components measured, n", "comps")
    r("Age, years, median (IQR)", "age")
    r("Gender code 1 of the released 0/1 variable, %", "gender1")
    r("Sensitised (database variable), %", "sens")
    r("Positive to ≥ 1 component, %", "pos1")
    r("Positive components per patient, median (IQR)", "breadth")
    r("Components positive per patient, %, median (IQR)", "pctpos")
    return rows


def summ():
    return json.load(open(OUT / "summary.json"))


def tableS2():
    """從晶片欄位到分析的成分對（流程表）。"""
    s = summ()
    bc = s["pairs_estimated_by_chip"]
    models = set(pd.read_csv(OUT / "component_plddt.csv").component)
    st = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    n_model, n_model_paired = len(models), len(models & (set(st.a) | set(st.b)))
    return [["Step", "n"],
            ["Allergen-specific IgE columns in the public file", f"{s['n_component_columns']}"],
            ["  Excluded: whole extracts or carbohydrate determinants", f"{s['n_extract_or_ccd']}"],
            ["  Excluded: listed in WHO/IUIS without a UniProt accession", f"{s['n_no_accession']}"],
            ["Single allergenic proteins mapped to UniProt", f"{s['n_mapped']}"],
            ["  With an AlphaFold model", f"{n_model}"],
            ["  Of these, in at least one estimable pair with another modelled component",
             f"{n_model_paired} ({n_model - n_model_paired} had fewer than 15 positive patients "
             f"on every platform)"],
            ["Component pairs with an adjusted co-sensitisation estimate, by platform",
             f"ISAC v1 {bc['ISAC_V1']:,}; ISAC v2 {bc['ISAC_V2']:,}; ALEX2 {bc['ALEX']:,}"],
            ["Pairs after pooling across platforms", f"{s['pairs_pooled']:,}"],
            ["  With an AlphaFold model for both members", f"{s['pairs_with_structure']:,}"],
            ["  Excluded: same source species (co-exposure)", f"{s['excl_same_species']:,}"],
            ["  Excluded: domains cut from one precursor (Hev b 6.01/6.02)", f"{s['excl_subchain']:,}"],
            ["  Excluded: a member shorter than 50 residues", f"{s['excl_short']:,}"],
            ["  Excluded: a member with no Pfam or SCOP superfamily annotation", f"{s['excl_unannotated']:,}"],
            ["Pairs analysed", f"{s['pairs_core']:,}"],
            ["  Homologous (shared Pfam family or SCOP superfamily)", f"{s['pairs_homologous']:,}"],
            ["  Non-homologous", f"{s['pairs_nonhomologous']:,}"],
            ["Components represented in the analysed pairs", f"{s['components_in_core']}"]]


def _disp(x):
    from common import display_name
    return display_name(x)


def _band_limits(band):
    nums = [float(v) for v in re.findall(r"\d\.\d+", band)]
    if band.startswith("<"):
        return -1, nums[0]
    if band.startswith("\u2265"):
        return nums[0], 2
    return nums[0], nums[1]


def table2_main():
    """同源對：TM-score 區間 → 共同致敏強度（臨床可查的對照表）。

    報 median OR 與 95% CI（節點自助法）；IQR 描述離散度，推論要用 CI。
    代表配對取「log OR 最接近該帶中位數」者——規則決定，不是手挑。
    """
    b = json.load(open(OUT / "bands.json"))
    t = pd.read_csv(OUT / "tm_bands_ci.csv")
    rows = [["TM-score", "Pairs, n", "Median odds ratio (95% CI)", "Pairs with z > 3, %",
             "Pair at the median of the band (TM-score)"]]
    import figs
    c = figs.core()
    tm = {}
    for x in c.itertuples():
        tm[(_disp(x.a), _disp(x.b))] = tm[(_disp(x.b), _disp(x.a))] = x.tm
    for r in t.itertuples():
        a_, b_ = r.example.split(" \u2013 ")
        v = tm[(a_, b_)]
        lo_, hi_ = _band_limits(r.band)
        assert lo_ <= v < hi_, (r.band, r.example, v)      # 代表配對一定要落在自己的帶內
        rows.append([r.band, f"{r.n}",
                     f"{r.median_or:.1f} ({r.lo:.1f} to {r.hi:.1f})",
                     f"{r.pct_z3:.0f}", f"{r.example} ({v:.3f})"])
    p_trend = b["trend_p_homologous"]
    rows.append(["Trend across bands", "", f"P {'< 0.001' if p_trend < 0.001 else f'= {p_trend:.3f}'}",
                 "", ""])
    return rows


def tableS3():
    """判別力與重分類。"""
    dd = json.load(open(OUT / "discrimination.json"))
    rows = [["Measure", "Value (95% CI)"]]
    rows.append(["Pairs analysed", f"{dd['n_core']:,}"])
    rows.append([f"Strongly co-sensitised (z > {dd['z_threshold']:.0f})",
                 f"{dd['n_strong']:,} ({100 * dd['n_strong'] / dd['n_core']:.1f}%)"])
    cd = json.load(open(OUT / "coef_diff.json"))
    rows.append(["  of which not homologous",
                 f"{cd['strong_nonhom']:,} ({100 * cd['strong_nonhom'] / cd['strong']:.0f}%)"])
    rows.append(["Area under the ROC curve", ""])
    names = {"tm": "  TM-score", "seqid80": "  Best 80-aa window identity",
             "seqid": "  Full-length sequence identity"}
    for r in dd["auc_all"]:
        rows.append([names[r["predictor"]], f"{r['auc']:.2f} ({r['lo']:.2f} to {r['hi']:.2f})"])
    for r in dd["auc_homologous"]:
        rows.append([names[r["predictor"]] + ", homologous pairs only",
                     f"{r['auc']:.2f} ({r['lo']:.2f} to {r['hi']:.2f})"])
    rows.append(["Binary criteria", ""])
    lab = {"codex": "  Sequence criterion (\u2265 35% over 80 residues)",
           "struct": "  Structural threshold (TM-score \u2265 0.5)"}
    for r in dd["points_all"]:
        rows.append([lab[r["criterion"]] + ": pairs flagged", f"{r['flagged']:,}"])
        rows.append(["    Sensitivity",
                     f"{r['sens']:.2f} ({r['sens_lo']:.2f} to {r['sens_hi']:.2f})"])
        rows.append(["    Specificity",
                     f"{r['spec']:.3f} ({r['spec_lo']:.3f} to {r['spec_hi']:.3f})"])
    rows.append(["Pairs by the two criteria", ""])
    nm = {(True, True): "  Both criteria met", (True, False): "  Structural threshold only",
          (False, True): "  Sequence criterion only", (False, False): "  Neither"}
    order = [(True, True), (True, False), (False, True), (False, False)]
    cells = {(r["structural"], r["codex"]): r for r in dd["reclassification"]}
    for k in order:
        r = cells[k]
        rows.append([nm[k], f"{r['n']:,} pairs; {r['strong']:,} strongly co-sensitised "
                            f"({r['pct_strong']:.0f}%); median odds ratio {r['median_or']:.1f}"])
    return rows


def tableS4():
    """結構達標、序列不達標的那 116 對，全部列出。"""
    q = pd.read_csv(OUT / "quadrant_struct_only.csv")
    rows = [["Allergen pair", "Source species", "TM-score", "80-aa identity, %",
             "Odds ratio (95% CI)", "z"]]
    for r in q.itertuples():
        lo, hi = np.exp(r.beta - 1.96 * r.se), np.exp(r.beta + 1.96 * r.se)
        rows.append([f"{disp(r.a)} \u2013 {disp(r.b)}",
                     f"<i>{r.org_a}</i> / <i>{r.org_b}</i>",
                     f"{r.tm:.2f}", f"{r.seqid80:.0f}",
                     _or_ci(np.exp(r.beta), lo, hi), f"{r.z:.1f}"])
    return rows


def tableS5():
    """對半分樣本：結局改在另一半病人身上定義。"""
    h = json.load(open(OUT / "holdout.json"))
    f = lambda k, nd=2: (f"{h[k]['median']:.{nd}f} ({h[k]['lo']:.{nd}f} to "
                         f"{h[k]['hi']:.{nd}f})")
    rows = [["Quantity", "Median across splits (range)"]]
    rows.append(["Splits", f"{h['reps']}"])
    rows.append(["Pairs estimable in both halves", f"{h['n_pairs']:,}"])
    rows.append(["  of which homologous", f"{h['n_homologous']}"])
    rows.append(["Agreement between the two halves", ""])
    rows.append(["  Spearman \u03c1 between the two halves' log odds ratios",
                 f('rho_between_halves')])
    sa = h["sign_agreement"]
    rows.append(["  Pairs with the same sign in both halves, %",
                 f"{100 * sa['median']:.0f}% ({100 * sa['lo']:.0f} to {100 * sa['hi']:.0f})"])
    rows.append(["Gradient within homologous pairs", ""])
    rows.append(["  Spearman \u03c1 between TM-score and co-sensitisation, each half alone",
                 f('rho_tm')])
    rows.append(["Discrimination with the outcome defined in the other half", ""])
    rows.append(["  Area under the curve, TM-score", f('auc_tm_out_of_sample')])
    rows.append(["  Area under the curve, 80-aa window identity",
                 f('auc_seq_out_of_sample')])
    rows.append(["Pairs above the structural threshold and below the sequence criterion", ""])
    rows.append(["  Pairs", f('struct_only_n', 0)])
    rows.append(["  Strongly co-sensitised in the other half", f('struct_only_strong', 0)])
    npct = h["neither_pct"]
    rows.append(["  Strongly co-sensitised among pairs meeting neither criterion, %",
                 f"{npct['median']:.1f}% ({npct['lo']:.1f} to {npct['hi']:.1f})"])
    return rows


def tableS6():
    """審閱要求的補強敏感度分析（analysis/extra.py）。"""
    e = json.load(open(OUT / "extra.json"))
    f2 = lambda v: f"{v:.2f}"
    rows = [["Analysis", "Result"]]
    rows.append(["Definition of strong co-sensitisation", ""])
    for r in e["z_thresholds"]:
        rows.append([f"  z > {r['z']:.0f}: pairs strongly co-sensitised",
                     f"{r['n_strong']:,} ({r['pct']:.1f}%); above the structural threshold and "
                     f"below the sequence criterion {r['struct_only_strong']} of "
                     f"{r['struct_only_n']}, against {r['neither_pct']:.1f}% of pairs meeting "
                     f"neither"])
    rows.append(["Multiplicity", ""])
    fdr = e["fdr"]
    rows.append(["  Pooled estimates with a Benjamini-Hochberg q below 0.05",
                 f"{fdr['q05']:,} of {fdr['n']:,}; all {fdr['strong_q05']:,} strongly "
                 f"co-sensitised pairs qualify"])
    rows.append(["  Of the pairs flagged by structure alone",
                 f"{fdr['struct_only_q05']} of {fdr['struct_only_n']} have q < 0.05"])
    rows.append(["Comparison of the two measures", ""])
    da, db = e["delong_all"], e["auc_diff_boot_all"]
    rows.append(["  Difference in area under the curve, all pairs (structure minus sequence)",
                 f"{da['diff']:.3f}; DeLong P = {da['p']:.2f}; bootstrap over allergens "
                 f"{db['lo']:.2f} to {db['hi']:.2f}"])
    dh, dbh = e["delong_homologous"], e["auc_diff_boot_homologous"]
    rows.append(["  Difference in area under the curve, homologous pairs",
                 f"{dh['diff']:.3f}; DeLong P = {dh['p']:.3f}; bootstrap over allergens "
                 f"{dbh['lo']:.2f} to {dbh['hi']:.2f}"])
    rows.append(["Normalisation of the TM-score", ""])
    n = e["normalisation"]
    rows.append(["  Spearman rho within homologous pairs, longer chain (used throughout) "
                 "and shorter chain",
                 f"{n['rho_long']:.2f} and {n['rho_short']:.2f}"])
    rows.append(["  Non-homologous pairs reaching a TM-score of 0.5, longer and shorter chain",
                 f"{0} and {n['nonhom_ge05_short']} (highest TM-score "
                 f"{n['nonhom_max_long']:.2f} and {n['nonhom_max_short']:.2f})"])
    rows.append(["  Pairs above the structural threshold and below the sequence criterion",
                 f"{n['struct_only_long']} of which {n['struct_only_long_strong']} strongly "
                 f"co-sensitised, against {n['struct_only_short']} of which "
                 f"{n['struct_only_short_strong']} with the shorter chain"])
    rows.append(["  Pairs meeting the sequence criterion only",
                 f"{n['seq_only_long']} and {n['seq_only_short']}"])
    rows.append(["Continuous rather than dichotomised IgE", ""])
    c = e["continuous"]
    rows.append(["  Partial Spearman between the two components' IgE values given breadth, "
                 "against TM-score, homologous pairs",
                 f2(c["rho_tm_vs_continuous_homologous"])])
    rows.append(["  The same against 80-aa window identity",
                 f2(c["rho_seq_vs_continuous_homologous"])])
    rows.append(["  Agreement with the dichotomised estimate, all analysed pairs",
                 f2(c["agreement_with_binary"])])
    rows.append(["Breadth computed only from components of other families", ""])
    dj = e["disjoint_breadth"]
    rows.append(["  Spearman rho between TM-score and co-sensitisation, homologous pairs",
                 f2(dj["rho_tm_homologous"])])
    rows.append(["  Median odds ratio of non-homologous pairs",
                 f"{dj['median_or_nonhom']:.1f}"])
    rows.append(["  Agreement with the main estimate, all analysed pairs",
                 f2(dj["agreement_with_main"])])
    rows.append(["Heterogeneity between platforms", ""])
    h = e["heterogeneity"]
    rows.append([f"  Pairs estimated on more than one platform",
                 f"{h['n_pairs_multi']:,}; median I2 {h['median_I2']:.0f}%, "
                 f"{h['pct_I2_gt50']:.0f}% above 50%"])
    return rows


def tableS2_sens():
    """敏感度分析（ρ 一律附節點自助法 95% CI）。"""
    s = summ()
    mt = {(r["stratum"]): r for r in s["mantel"]}
    ci = json.load(open(OUT / "rho_ci.json"))
    pl = json.load(open(OUT / "plddt.json"))
    rb = s["robust"]
    get = lambda an, t: next(r for r in rb if r["analysis"] == an and r["term"] == t)
    f = lambda k: f"{ci[k]['rho']:.2f} ({ci[k]['lo']:.2f} to {ci[k]['hi']:.2f})"
    f3 = lambda k: (f"{ci[k]['rho']:.3f} ({ci[k]['lo']:.3f} to "
                    f"{ci[k]['hi']:.3f})").replace("-", "\u2212")
    # 欄名只寫預測變數；量是什麼寫在分段列——聯合模型的係數不是 ρ，不可放在「Spearman ρ」欄名下
    rows = [["Analysis", "Pairs, n", "TM-score", "80-aa window identity"],
            ["Spearman \u03c1 with co-sensitisation (95% CI for TM-score)", "", "", ""]]
    hs = mt["homologous, sequence identity"]
    rows.append(["Homologous pairs, main analysis", f"{ci['homologous']['n']}",
                 f("homologous"), f"{hs['spearman']:.2f}"])
    fl = mt["homologous, full-length identity"]
    rows.append(["  Full-length identity instead of the 80-aa window", f"{ci['homologous']['n']}",
                 "", f"{fl['spearman']:.2f}"])
    xx = get("collinearity", "Spearman TM vs identity")
    rows.append(["  TM-score against 80-aa identity (correlation between the two measures)",
                 f"{int(xx['n'])}", f"{xx['estimate']:.2f}", ""])
    for key, lab in [("isac_v1", "ISAC v1"), ("isac_v2", "ISAC v2"), ("alex", "ALEX2")]:
        q = get("by chip", {"isac_v1": "ISAC_V1", "isac_v2": "ISAC_V2",
                            "alex": "ALEX"}[key] + " rho identity")
        rows.append([f"  Estimated on {lab} only", f"{ci[key]['n']}", f(key),
                     f"{q['estimate']:.2f}"])
    q = get("different taxonomic order", "rho identity")
    rows.append(["  Source species from different taxonomic orders",
                 f"{ci['different_order']['n']}", f("different_order"),
                 f"{q['estimate']:.2f}"])
    rows.append([f"  Both AlphaFold models with mean pLDDT \u2265 {pl['cut']:.0f}",
                 f"{ci['plddt70']['n']}", f("plddt70"), ""])
    lo = get("leave one homology group out", "rho TM min")
    hi = get("leave one homology group out", "rho TM max")
    rows.append([f"  Leaving out each of {int(lo['n'])} homology groups in turn", "",
                 f"{lo['estimate']:.2f} to {hi['estimate']:.2f}", ""])
    cx = mt["homologous, below Codex sequence criterion"]
    rows.append(["  Pairs below the sequence criterion (< 35% identity over 80 residues)",
                 f"{ci['below_codex']['n']}",
                 f("below_codex") + f", P = {cx['p_perm']:.3f}", ""])
    n = mt["non-homologous"]
    rows.append(["Non-homologous pairs, main analysis", f"{ci['nonhomologous']['n']:,}",
                 f3("nonhomologous") + f", P = {n['p_perm']:.2f}", ""])
    fr = mt["non-homologous, excluding fragments"]
    rows.append(["  Excluding fragment sequences", f"{ci['nonhomologous_nofrag']['n']:,}",
                 f3("nonhomologous_nofrag") + f", P = {fr['p_perm']:.2f}", ""])
    cu = mt["non-homologous, hand-curated family labels"]
    rows.append(["  Homology defined by curated WHO/IUIS family names instead",
                 f"{ci['nonhomologous_curated']['n']:,}",
                 f3("nonhomologous_curated") + f", P = {cu['p_perm']:.2f}", ""])
    jt = get("joint model", "TM (rank, SD)")
    js = get("joint model", "sequence identity (rank, SD)")
    rows.append(["Both predictors in one model, homologous pairs: regression coefficient per SD "
                 "of rank (not \u03c1)", "", "", ""])
    cd = json.load(open(OUT / "coef_diff.json"))
    ci_ = lambda v, lo_, hi_: (f"{v:.2f} ({lo_:.2f} to {hi_:.2f})").replace("-", "\u2212")
    # 聯合模型報節點自助法 CI，不報置換 P：同源對很稀疏，計入的置換中位數只有 35 對可比，
    # 虛無分布極寬，係數的置換 P 沒有檢定力（曾印出 P = 0.500，與 CI 0.08–1.07 矛盾）
    rows.append(["  Coefficient (95% CI)", f"{cd['n']}", ci_(cd["b_tm"], *cd["b_tm_ci"]),
                 ci_(cd["b_seq"], *cd["b_seq_ci"])])
    rows.append(["  Difference, structure minus sequence (95% CI)",
                 f"{cd['n']}", f"{cd['diff']:.2f} ({cd['diff_ci'][0]:.2f} to "
                 f"{cd['diff_ci'][1]:.2f})".replace("-", "\u2212"), ""])
    return rows


def items():
    s = summ()
    import figs
    p_all = figs._pairs()[0]
    n_sub = int((p_all.same_species & p_all.subchain).sum())
    cp = pd.read_csv(OUT / "component_plddt.csv").set_index("component").plddt
    hh = figs.core(); hh = hh[hh.homologous]
    lowm = hh[(hh.a.map(cp) < 70) | (hh.b.map(cp) < 70)]
    who = pd.concat([lowm.a[lowm.a.map(cp) < 70], lowm.b[lowm.b.map(cp) < 70]]).value_counts()
    s2 = dict(n_below=len(lowm), top=_disp(who.index[0]), n_top=int(who.iloc[0]),
              top_plddt=cp[who.index[0]])
    dd = json.load(open(OUT / "discrimination.json"))
    pl = json.load(open(OUT / "plddt.json"))
    q = dd["quadrant_struct_only"]
    bg = dd["quadrant_neither"]
    main = [
        dict(kind="table",
             caption="Table 1. Patients contributing allergen-chip results, by platform.",
             rows=table1(),
             note="Allergen Chip Challenge public release (11 French university hospitals, "
                  "2014-2023). A component counts as measured on a platform when a value is "
                  "present for at least 80% of that platform's patients; positivity is 0.3 ISU-E "
                  "or kUA/L. The last column counts components measured on at least one platform. "
                  "Sex is released as an undocumented 0/1 code, so the two levels cannot be "
                  "labelled; the distribution is given as released, and 62 patients have no "
                  "value. "
                  "IQR, interquartile range."),
        dict(kind="table",
             caption="Table 2. Co-sensitisation between homologous allergens from different "
                     "species, by predicted structural similarity.",
             rows=table2_main(),
             note=f"Homologous pairs share a Pfam family or a SCOP superfamily "
                  f"(n = {s['pairs_homologous']}). TM-score compares AlphaFold models of the two "
                  "proteins (1 = identical fold; values above 0.5 usually indicate the same fold). "
                  "Odds ratios are adjusted for the number of other components to which the "
                  "patient is positive and pooled across the three platforms; confidence "
                  "intervals are from 2,000 bootstrap resamples of allergens, not of pairs, "
                  "because pairs sharing a component are not independent. The pair shown for each "
                  "band is the one whose odds ratio lies closest to the median of that band. "
                  "The trend P value is from permutation of allergen labels. CI, confidence "
                  "interval.",),
        dict(kind="figure", path=OUT / "Figure_1.png",
             caption="Figure 1. Study design, derivation of the analysed allergen pairs, and the "
                     "effect of adjusting for sensitisation breadth.",
             note=f"(A) Analysis steps. (B) Derivation of the analysed pairs, with bar length proportional to the number of "
                  f"pairs; counts in grey are "
                  f"the pairs removed at each step (Supplementary Table S1). (C) Percentage of a "
                  f"platform's components to which each patient is positive; the last bin includes "
                  f"the few patients above 60%. (D) Unadjusted "
                  f"against breadth-adjusted log odds ratios for all {s['pairs_pooled']:,} "
                  f"component pairs; the dashed line is equality. Without adjustment almost every "
                  f"pair is positively associated (median log odds ratio "
                  f"{s['median_crude_logOR']:.2f}); after adjustment the median is "
                  f"{s['median_adjusted_logOR']:.2f}. OR, odds ratio."),
        dict(kind="figure", path=OUT / "Figure_2.png",
             caption="Figure 2. Predicted structural similarity grades co-sensitisation within, "
                     "but not beyond, homologous allergen families.",
             note=f"(A) Co-sensitisation for pairs from the same source species "
                  f"(n = {s['same_species_n']}; the {s['same_species_n'] + n_sub} same-species pairs "
                  f"in Supplementary Table S1 include {n_sub} pairs in which one member is a "
                  f"domain cut from one precursor (all involve Hev b 6.01), "
                  f"which are not shown), homologous pairs from different species "
                  f"(n = {s['homologous_n']}) and non-homologous pairs "
                  f"(n = {s['nonhomologous_n']:,}); bars give the median and interquartile range "
                  f"(median odds ratios {np.exp(s['same_species_median']):.1f}, "
                  f"{np.exp(s['homologous_median']):.1f} and "
                  f"{np.exp(s['nonhomologous_median']):.2f}). (B) Homologous pairs against "
                  f"TM-score; filled points share at least 35% identity over an 80-residue "
                  f"window, the Codex Alimentarius screening criterion, and point area is "
                  f"proportional to 1/SE of the estimate; red segments are the medians "
                  f"of the bands in Table 2. (C) Non-homologous pairs, shaded by density; the "
                  f"highest TM-score reached was {s['nonhomologous_tm_max']:.2f}, so the axis "
                  f"stops at 0.56; the vertical axis is kept the same as in A and B so that the "
                  f"three panels can be compared directly. (D) Median odds ratio of each band with its 95% confidence "
                  f"interval. \u03c1, Spearman correlation; P values are from permutation of "
                  f"allergen labels."),
        dict(kind="figure", path=OUT / "Figure_3.png",
             caption="Figure 3. What the sequence screening criterion and a structural threshold "
                     "each identify.",
             note=f"All {dd['n_core']:,} analysed pairs. Strong co-sensitisation is z > 3 "
                  f"({dd['n_strong']:,} pairs, {100 * dd['n_strong'] / dd['n_core']:.0f}%). "
                  f"(A) Each pair by sequence and structural similarity, coloured by its adjusted "
                  f"log odds ratio; dashed lines are the Codex Alimentarius criterion (35% "
                  f"identity over 80 residues) and a TM-score of 0.5, each quadrant is labelled "
                  f"with its number of pairs and median odds ratio, and pairs with no detectable "
                  f"local sequence similarity are plotted at \"none\". Open circles are the pairs "
                  f"shown in C; two of the 14 have no detectable sequence similarity and nearly "
                  f"the same TM-score, so they overlap at \"none\". (B) Pairs and strongly "
                  f"co-sensitised pairs in each cell of the "
                  f"two criteria. (C) The 14 pairs with the highest z among the {q['n']} that "
                  f"fall below the sequence criterion but reach a TM-score of 0.5, with their "
                  f"TM-score and 80-aa window identity; all {q['n']} are listed in Supplementary "
                  f"Table S3. (D) Discrimination for strong co-sensitisation, with the two binary "
                  f"criteria as points; confidence intervals are from 2,000 bootstrap resamples "
                  f"of allergens. The two similarity measures discriminate equally well; in a joint model on "
                  f"homologous pairs the two coefficients cannot be distinguished (Supplementary "
                  f"Table S4). The evidence that structure adds "
                  f"information is the {q['n']} pairs it flags on its own (B, C) and their "
                  f"replication in the split-sample analysis (Supplementary Table S6), not a "
                  f"difference in the areas under the curve or in the joint-model coefficients. AUC, area under the "
                  f"receiver operating characteristic curve; CI, confidence interval."),
        dict(kind="figure", path=OUT / "Figure_4.png",
             caption="Figure 4. Strong co-sensitisation reproduces the recognised allergen "
                     "families, and a shared fold can carry it where sequence cannot.",
             note=f"(A) Components joined when their pair is strongly co-sensitised (z > 3), "
                  f"placed by a force-directed layout seeded by family; nodes are coloured for "
                  f"the eight families with the most components and edges between homologous "
                  f"components are drawn darker. No structural information was used to draw this "
                  f"panel. (B) Median adjusted log odds ratio for each pair of families, where at "
                  f"least three pairs contribute. The colour scale is centred at zero and is "
                  f"stretched differently on either side (\u22122 to 0 and 0 to 5); values above 5, "
                  f"such as the profilin and tropomyosin diagonal cells, take the end colour, and "
                  f"every cell is labelled with its value. The 7S and 11S seed storage globulins are "
                  f"positive off the diagonal as well as on it, which is expected because they "
                  f"belong to one cupin superfamily. (C, D) AlphaFold models superposed by "
                  f"TM-align, drawn as cartoons (first allergen dark blue, second amber); signal "
                  f"peptides annotated in UniProt, residues with low model confidence "
                  f"(pLDDT < 70) and fragments shorter than five residues left by that removal "
                  f"are not drawn; the models are trimmed for display only, and TM-scores "
                  f"are computed on the complete models. Cry j 1, a pectate lyase of Japanese cedar, and "
                  f"Pla a 2, a polygalacturonase of London plane, belong to different Pfam "
                  f"families but share the pectin lyase-like superfamily (a right-handed "
                  f"\u03b2-helix); they have no detectable sequence similarity and the same fold; Ara h 2 and Api g 2 belong to one superfamily but their folds "
                  f"have diverged. The TM-score shown is the one used throughout, normalised by "
                  f"the longer chain."),
    ]
    supp = [
        dict(kind="table", caption="Supplementary Table S1. From chip columns to analysed "
                                   "allergen pairs.",
             rows=tableS2(),
             note="Allergen names were resolved against the WHO/IUIS Allergen Nomenclature "
                  "database, which lists a UniProt accession for each isoallergen. A pair was "
                  "The exclusions are applied in the order listed and each pair is counted once "
                  "at the first rule it meets, so the five same-species pairs that also involve a "
                  "domain of one precursor are counted as same-species. A pair was "
                  "estimated on a platform when each member had at least 15 positive patients, at "
                  "least 50 patients had results for both and at least 5 were positive to both. "
                  "Counts of components are chip analytes: a protein measured as more than one "
                  "analyte counts once per analyte (Ole e 7 is measured on ISAC and, as a "
                  "research-use analyte marked RUO, on ALEX2)."),
        dict(kind="table", caption="Supplementary Table S2. Discrimination for strong "
                                   "co-sensitisation and pairs flagged by each criterion.",
             rows=tableS3(),
             note=f"Strong co-sensitisation is a z above {dd['z_threshold']:.0f} for "
                  "the pooled adjusted estimate. These are measures of agreement within the same "
                  "cohort, not external validation: the outcome is derived from the same odds "
                  "ratios. Confidence intervals are from 2,000 bootstrap resamples of allergens, "
                  "not of pairs. Full-length identity is the number of identical residues in the "
                  "best local alignment divided by the length of the longer sequence. ROC, "
                  "receiver operating characteristic."),
        dict(kind="table", caption="Supplementary Table S3. All pairs that fall below the "
                                   "sequence criterion but reach a TM-score of 0.5.",
             rows=tableS4(),
             note=f"Ordered by z. {q['strong']} of the {q['n']} pairs are strongly co-sensitised "
                  f"(z > 3) and the median odds ratio is {q['median_or']:.1f}, against "
                  f"{100 * bg['strong'] / bg['n']:.0f}% and {bg['median_or']:.1f} among the "
                  f"{bg['n']:,} pairs that meet neither criterion. CI, confidence interval."),
        dict(kind="table", caption="Supplementary Table S4. Sensitivity analyses.",
             rows=tableS2_sens(),
             note="P values are from 10,000 permutations of allergen labels, because pairs that "
                  "share a component are not independent; a permutation counts only when at least "
                  "30 of the analysed pairs remain comparable. Among homologous pairs a permutation "
                  "leaves a median of 35 comparable pairs, so these P values are conservative, "
                  "and the joint model is reported with confidence intervals only (5,000 bootstrap "
                  "resamples of allergens). Other confidence intervals are from 2,000 bootstrap "
                  "resamples of allergens. Full-length identity is the number of identical residues "
                  "in the best local alignment divided by the length of the longer sequence. "
                  "pLDDT is the AlphaFold per-residue confidence, averaged "
                  f"over the model; {pl['n_components_below']} of {pl['n_components']} components "
                  "have a mean pLDDT below 70."),
        dict(kind="table", caption="Supplementary Table S5. Additional sensitivity analyses.",
             rows=tableS6(),
             note="Benjamini-Hochberg q values are computed over all 18,848 pooled estimates. "
                  "DeLong's test treats pairs as independent observations, which they are not, "
                  "so the bootstrap over allergens is the comparison we rely on. The TM-score is "
                  "normalised by the longer chain throughout the article; normalising by the "
                  "shorter chain lets a short protein resemble a much longer one, and the "
                  "non-homologous pairs it brings above 0.5 have a median shorter-chain length of "
                  "85 residues and are not co-sensitised. Breadth computed from components of "
                  "other families removes the arithmetic dependence created by conditioning on a "
                  "count that includes the two components being compared. I2 is the proportion of "
                  "variance between platform-specific estimates not explained by sampling error."),
        dict(kind="table", caption="Supplementary Table S6. Split-sample analysis with the "
                                   "outcome defined in a disjoint set of patients.",
             rows=tableS5(),
             note="Patients were split in half at random within each platform and the "
                  "co-sensitisation of every pair was re-estimated independently in each half, "
                  "repeated five times. Because TM-score and sequence identity use no patient "
                  "data, defining strong co-sensitisation in one half and evaluating the "
                  "similarity measures against it removes the shared estimation that makes the "
                  "main discrimination analysis internal. Fewer pairs are estimable than in the "
                  "main analysis because each half has half the patients; the areas under the "
                  "curve are higher than in the full cohort (Supplementary Table S2) because "
                  f"the estimable pairs are a different and smaller set "
                  f"({json.load(open(OUT / 'holdout.json'))['n_pairs']:,} rather than "
                  f"{dd['n_core']:,})."),
        dict(kind="figure", path=OUT / "Figure_S1.png",
             caption="Supplementary Figure S1. Sensitivity analyses.",
             note="(A) Spearman correlation between TM-score and adjusted co-sensitisation, with "
                  "95% confidence intervals from 2,000 bootstrap resamples of allergens. Dark "
                  "points are the two main analyses; the rows beneath each repeat it under a "
                  "different restriction (Supplementary Table S4). (B) The same correlation among "
                  "homologous pairs after removing each homology group in turn; the red line is "
                  "the value for all homologous pairs."),
        dict(kind="figure", path=OUT / "Figure_S2.png",
             caption="Supplementary Figure S2. AlphaFold model confidence.",
             note=f"(A) Mean pLDDT of the model of each component; the line marks 70, the usual "
                  f"threshold for a reliable backbone. (B) TM-score of each homologous pair "
                  f"against the lower of the two models' mean pLDDT. Models with lower confidence "
                  f"do reach lower TM-scores (\u03c1 = {pl['rho_plddt_tm']:.2f}), but restricting "
                  f"to pairs in which both models reach 70 leaves the gradient unchanged "
                  f"({pl['rho_all']:.2f} to {pl['rho_high_confidence']:.2f}). Of the "
                  f"{s2['n_below']} homologous pairs with a model below 70, {s2['n_top']} involve "
                  f"{s2['top']}, whose model has a mean pLDDT of {s2['top_plddt']:.1f}, so they lie "
                  f"against the line; the line is drawn behind the points."),
        dict(kind="figure", path=OUT / "Figure_S3.png",
             caption="Supplementary Figure S3. Null distributions of the permutation test.",
             note="Spearman correlation recomputed after permuting the allergen labels, which "
                  "preserves the dependence between pairs that share a component; 10,000 "
                  "permutations each. A permutation counts only when at least 30 of the pairs "
                  "remain comparable after relabelling. Among homologous pairs a median of 35 do, "
                  "which is why the null distribution in A is wide and the test conservative. The "
                  "observed value is marked in red."),
        dict(kind="figure", path=OUT / "Figure_S4.png",
             caption="Supplementary Figure S4. The gradient estimated separately on each "
                     "platform.",
             note="Adjusted log odds ratio from each platform alone against TM-score; homologous "
                  "pairs are coloured, non-homologous pairs are grey. Pairs are fewer on each "
                  "platform than in the pooled analysis because the platforms measure different "
                  "components."),
    ]
    return main, supp


# ------------------------------------------------------------------ source data
def _df(rows):
    return pd.DataFrame(rows[1:], columns=rows[0])


def source_data(path, main, supp):
    """全部圖表一份 Excel：表格照 Word 內容，圖每個面板一張分頁（source_data.py）。"""
    import source_data as SD
    figs_ = SD.panels()
    tabs = {"Table 1": main[0], "Table 2": main[1]}
    tabs.update({f"Table S{i}": it for i, it in enumerate(
        [x for x in supp if x["kind"] == "table"], 1)})
    t = lambda k: (k, tabs[k]["caption"], _df(tabs[k]["rows"]))
    order = [t("Table 1"), t("Table 2")]
    for k in "1234":
        order += [x for x in figs_ if x[0].split()[1][0] == k]
    order += [t(f"Table S{i}") for i in range(1, 7)]
    order += [x for x in figs_ if x[0].split()[1].startswith("S")]
    assert len(order) == len(figs_) + 8, "有分頁沒排進去"
    SD.write(path, order)


def data_s1(path):
    """補充資料 S1：全部成分對的校正後共同致敏估計＋結構／序列相似度（本文的可重用資源）。

    附 UniProt accession 與 Pfam／SCOP 註解，讀者才能自行核對同源判定。
    """
    import homology
    adj = pd.read_csv(OUT / "cosens_pooled.csv")
    import structure_cosens as S
    st = pd.read_csv(OUT / "pairs_structure_cosens.csv")[["a", "b", "tm", "seqid", "seqid80",
                                                          "same_species", "homologous", "subchain",
                                                          "len_a", "len_b", "annotated"]]
    d = adj.merge(st, on=["a", "b"], how="left")
    # 每一對為什麼進或不進主分析——讀者要能從這個檔重現 8,239 對（順序同 Table S1）
    reason = np.select(
        [d.tm.isna(), d.same_species == True, d.subchain == True,
         (d.len_a < S.MIN_LEN) | (d.len_b < S.MIN_LEN), d.annotated != True],
        ["no AlphaFold model for one member", "same source species",
         "domains cut from one precursor", "member shorter than 50 residues",
         "no Pfam or SCOP superfamily annotation"], default="")
    d["analysed"] = reason == ""
    d["exclusion_reason"] = reason
    assert d.analysed.sum() == len(figs_core := __import__("figs").core())
    assert (d.analysed & (d.homologous == True)).sum() == int(figs_core.homologous.sum())
    d = d.drop(columns=["subchain", "len_a", "len_b", "annotated"])
    m = families.table()
    m = m[(m.status == "ok") & m.accession.notna()]
    acc = m.groupby("component").accession.apply(lambda v: ";".join(sorted(set(v)))).to_dict()
    ann = homology.table()
    fam = m.drop_duplicates("component").set_index("component")
    q = pd.read_csv(OUT / "qvalues.csv")
    d = d.merge(q, on=["a", "b"], how="left")
    d["odds_ratio"] = np.exp(d.beta)
    d["ci_low"] = np.exp(d.beta - 1.96 * d.se)
    d["ci_high"] = np.exp(d.beta + 1.96 * d.se)
    for side in "ab":
        d[f"uniprot_{side}"] = d[side].map(acc)
        d[f"family_by_IUIS_name_{side}"] = d[side].map(fam.family)
        d[f"organism_{side}"] = d[side].map(fam.organism)
        d[f"pfam_scop_{side}"] = d[side].map(lambda c: ";".join(sorted(ann.get(c, []))))
    for side in "ab":
        d[side] = d[side].map(disp)      # `_RUO` 是平台標記，不是過敏原名稱的一部分
    d = d.rename(columns={"a": "component_a", "b": "component_b", "beta": "adjusted_log_OR",
                          "se": "SE", "k_chip": "platforms", "n_total": "patients",
                          "tm": "TM_score", "seqid": "full_length_identity_pct",
                          "seqid80": "best_80aa_window_identity_pct",
                          "p_value": "p_value_two_sided",
                          "q_value": "q_value_benjamini_hochberg"})
    d.to_csv(path, index=False)
    return len(d)


def main():
    dest = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DEST
    os.makedirs(dest, exist_ok=True)
    main_items, supp = items()
    for it in main_items + supp:
        if it["kind"] == "figure" and not os.path.exists(it["path"]):
            raise SystemExit(f"缺圖：{it['path']}（先跑 figs.py）")
    import meta
    out = os.path.join(dest, "Tables_and_Figures.docx")
    print(build(main_items, out, supp))
    meta.clean_docx(out)
    old = os.path.join(dest, "Supplementary_Figures_and_Tables.docx")
    if os.path.exists(old):
        os.remove(old)     # 已併入同一份 Word
    for it in main_items + supp:
        if it["kind"] == "figure":
            shutil.copy(it["path"], os.path.join(dest, os.path.basename(it["path"])))
            print("  png ->", os.path.basename(it["path"]))
    source_data(os.path.join(dest, "Source_data.xlsx"), main_items, supp)
    meta.clean_xlsx(os.path.join(dest, "Source_data.xlsx"))
    print("  xlsx -> Source_data.xlsx")
    print("  csv  -> Supplementary_Data_S1.csv", data_s1(os.path.join(dest, "Supplementary_Data_S1.csv")))


if __name__ == "__main__":
    main()
