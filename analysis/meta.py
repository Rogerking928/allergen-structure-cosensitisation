"""交付檔的文件屬性：作者寫使用者本人，清掉產生工具留下的字樣（2026-10-02 使用者要求）。"""
import datetime

AUTHOR = "Yen-Hsiang Wang"


def clean_docx(path):
    from docx import Document
    d = Document(path)
    c = d.core_properties
    c.author = AUTHOR
    c.last_modified_by = AUTHOR
    c.comments = ""
    c.created = c.modified = datetime.datetime.now()
    d.save(path)


def clean_xlsx(path):
    import openpyxl
    w = openpyxl.load_workbook(path)
    w.properties.creator = AUTHOR
    w.properties.lastModifiedBy = AUTHOR
    w.save(path)
    # openpyxl 會把自己寫進 docProps/app.xml 的 Application，改寫成 Excel 的預設字樣
    import os, re, shutil, tempfile, zipfile
    tmp = tempfile.mktemp(suffix=".xlsx", dir=os.path.dirname(path))
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "docProps/app.xml":
                data = re.sub(rb"<Application>[^<]*</Application>", b"<Application>Microsoft Excel</Application>", data)
                data = re.sub(rb"<AppVersion>[^<]*</AppVersion>", b"<AppVersion>16.0300</AppVersion>", data)
            zout.writestr(item, data)
    shutil.move(tmp, path)
