"""Figure 4C/D 的 cartoon 疊合圖（PyMOL 光線追蹤），取代原本的 Cα 軌跡。

要用裝了 pymol-open-source 的 Python 跑（系統 Python 沒有 PyMOL）：
  <pymol env>/bin/python analysis/overlay_render.py

疊合矩陣沿用 overlay.py 存下來的 foldseek TM-align 結果（u, t：x' = u·x + t），
套用後重算 Cα 與 overlay.py 存的座標比對，對不上就停——避免畫出跟 RMSD 檢查不同的疊合。
只畫 UniProt 的成熟鏈（Chain 範圍）且 pLDDT ≥ 70 的殘基，訊號肽不畫：Ara h 2 的訊號肽在 AlphaFold 裡是一條
長螺旋，會佔掉半張圖，而 IgE 遇到的是成熟蛋白。
"""
import json
import pathlib
import numpy as np
from pymol import cmd

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "out"
AFDB = ROOT / "data" / "afdb"
UNI = ROOT / "data" / "uniprot"
PAIRS = [("Cry_j_1", "Pla_a_2"), ("Ara_h_2", "Api_g_2")]
NAVY, AMBER = "0x31485d", "0xb8791f"
PLDDT_MIN = 70


def chain_range(acc):
    d = json.load(open(UNI / f"{acc}.json"))
    ch = [f for f in d.get("features", []) if f["type"] == "Chain"]
    if not ch:
        return None
    return ch[0]["location"]["start"]["value"], ch[0]["location"]["end"]["value"]


def render(a, b, view, size=1800, png=None):
    d = np.load(OUT / f"overlay_{a}_{b}.npz")
    qa, ta = str(d["acc_q"]), str(d["acc_t"])
    cmd.reinitialize()
    cmd.load(str(AFDB / f"{qa}.pdb"), "mq")
    cmd.load(str(AFDB / f"{ta}.pdb"), "mt")
    u, t = d["u"], d["tr"]
    M = [*u[0], t[0], *u[1], t[1], *u[2], t[2], 0, 0, 0, 1]
    cmd.transform_selection("mt", M, homogenous=1)
    got = np.array(cmd.get_coords("mt and name CA"))
    err = np.abs(got - d["t"]).max()
    if err > 0.01:
        raise SystemExit(f"{a}-{b}: 轉換後座標與 overlay.py 差 {err:.2f} Å，停")
    for obj, acc in (("mq", qa), ("mt", ta)):
        r = chain_range(acc)
        if r:
            cmd.remove(f"{obj} and not resi {r[0]}-{r[1]}")
    # AlphaFold 的 B-factor 欄是 pLDDT；< 70 屬低信心（Ara h 2 的 56–97 號是無序迴圈），不畫
    cmd.remove(f"b < {PLDDT_MIN}")
    # 拿掉低信心區後殘下的零碎片段（< 5 個殘基）也不畫，否則會在圖邊浮著一小截
    for obj in ("mq", "mt"):
        res = sorted({int(a.resi) for a in cmd.get_model(f"{obj} and name CA").atom})
        runs, cur = [], [res[0]]
        for r in res[1:]:
            if r == cur[-1] + 1:
                cur.append(r)
            else:
                runs.append(cur); cur = [r]
        runs.append(cur)
        for run in runs:
            if len(run) < 5:
                cmd.remove(f"{obj} and resi {run[0]}-{run[-1]}")
    cmd.hide("everything")
    cmd.dss()
    cmd.show("cartoon")
    cmd.color(NAVY, "mq")
    cmd.color(AMBER, "mt")
    cmd.set("cartoon_transparency", 0.35, "mt")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_gap_cutoff", 0)   # 拿掉的殘基不要用虛線連起來
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.35)
    cmd.set("specular", 0.15)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_highlight_color", "grey70")
    cmd.orient("mq or mt")
    for ax, ang in view:
        cmd.turn(ax, ang)
    cmd.zoom("mq or mt", buffer=1.5)
    png = pathlib.Path(png or OUT / f"overlay_{a}_{b}.png")
    cmd.png(str(png), width=size, height=size, dpi=700, ray=1)
    print(f"{a}-{b} -> {png.name}（最大座標差 {err:.3f} Å）")


if __name__ == "__main__":
    render("Cry_j_1", "Pla_a_2", [])
    render("Ara_h_2", "Api_g_2", [])
