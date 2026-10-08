"""正文對圖表本的三條檢查（照 MRI 標示那篇的做法移植過來）。

1. **正文的每一個數字都要在圖表裡找得到**。讀者查不到的數字不准寫。
2. **圖表的編號要照正文第一次引用的順序**，不可跳號。
3. **每一個圖表都至少被引用一次**，主文與補充都算。

  python3 code/check_manuscript.py [code/manuscript.md]
"""
import os
import pathlib
import re
import sys
import zipfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEST = pathlib.Path(os.environ.get(
    "DELIVERABLE_DIR", "/mnt/c/Users/roger/Desktop/還沒動/過敏原結構_共同致敏"))
BOOK = DEST / "Tables_and_Figures.docx"

# 屬於稿件本身的體例、不是結果的數字
ALLOWED = {
    "2.1", "2.2", "2.3", "2.4", "3.1", "3.2", "3.3", "3.4", "3.5",
    "1", "2", "3", "4", "11031", "250", "2014", "2023", "2026", "95",
    "80", "35", "0.3", "50", "11", "3,500",
    # 引文裡別人的世代人數，不是本文的結果
    "3,000", "23,000",
}


def book_text():
    x = zipfile.ZipFile(BOOK).read("word/document.xml").decode()
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", re.sub(r"</w:p>", "\n", x)))


def body(md):
    text = md[md.index("\nIntroduction\n"):md.index("\nCompeting interests:")]
    text = re.sub(r"\[\d+(?:[,–-]\d+)*\]", "", text)
    return re.sub(r"\{[^}]*\}", "", text)


def main(src="code/manuscript.md"):
    md = (ROOT / src).read_text(encoding="utf-8")
    b = body(md)
    book = book_text()
    bad = []

    nums = {n.rstrip(".,") for n in
            re.findall(r"(?<![\w.,/])\d[\d,]*(?:\.\d+)?(?![\w])", b)}
    for n in sorted(nums - ALLOWED):
        if n not in book:
            bad.append(f"正文寫了 {n}，但任何圖表裡都沒有這個數字")

    PAT = {"Figure": r"(?<!Supplementary )Figure (\d)",
           "Table": r"(?<!Supplementary )Table (\d)",
           "Supplementary Figure": r"Supplementary Figures? S(\d)",
           "Supplementary Table": r"Supplementary Tables? S(\d)"}
    present = {"Figure": 4, "Table": 2, "Supplementary Figure": 4,
               "Supplementary Table": 5}
    for name, pat in PAT.items():
        seen, order = set(), []
        for m in re.finditer(pat, b):
            k = int(m.group(1))
            if k not in seen:
                seen.add(k)
                order.append(k)
        if order != sorted(order):
            bad.append(f"{name} 第一次引用的順序不對：{order}")
        missing = sorted(set(range(1, present[name] + 1)) - seen)
        if missing:
            bad.append(f"{name} {missing} 從頭到尾沒被引用")

    for fig, panels in {"Figure 1": "ABCD", "Figure 2": "ABCD",
                        "Figure 3": "ABCD", "Figure 4": "ABC"}.items():
        n = fig.split()[1]
        gone = [p for p in panels
                if not re.search(rf"{fig}{p}|Figure {n}{p}|{n}{p}[,;)\s]", b)]
        if gone:
            bad.append(f"{fig} 沒被引用到的分格：{', '.join(gone)}")

    words = len(re.sub(r"\s+", " ", b).split())
    for line in bad:
        print("  ", line)
    print(f"\n本文 {words:,} 字，檢查了 {len(nums)} 個數字")
    print("OK" if not bad else "上面要修")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
