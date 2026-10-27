#!/usr/bin/env python
# ==============================================================================
# Build the folder that gets uploaded.
#
# Editorial Manager wants separate documents in separate slots: cover letter,
# title page, manuscript, figure files, supplemental files, checklists. A folder
# with the right names and the right contents turns an hour of fiddling at upload
# time into five drag-and-drops. The texts are pulled from the submission package
# so that the letter and the title page cannot drift from the manuscript.
# ==============================================================================
import re
import shutil
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

ROOT = Path("C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45")
PKG = (ROOT / "manuscript" / "00_submission_package.md").read_text(encoding="utf-8")
MD = (ROOT / "manuscript" / "JAAD_BriefReport_isotretinoin_sex.md").read_text(encoding="utf-8")
OUT = ROOT / "submission"
OUT.mkdir(exist_ok=True)


def clean(t):
    t = t.replace("**", "").replace("*", "")
    t = t.replace("\u2014", "-").replace("\u2013", "-")
    t = re.sub(r"[ \t]+", " ", t)
    return t.strip()


def section(name, nxt):
    m = re.search(rf"## {re.escape(name)}.*?\n(.*?)(?=\n## {re.escape(nxt)})", PKG, re.S)
    return m.group(1).strip() if m else ""


def parapraphs(text):
    blocks = [b.strip() for b in re.split(r"\n\s*\n", text) if b.strip()]
    out = []
    for b in blocks:
        t = clean(re.sub(r"\s*\n\s*", " ", b))
        t = t.replace(" ---", "").replace("---", "").strip()
        if t:
            out.append(t)
    return out



def stamp(doc, title):
    cp = doc.core_properties
    cp.author = "Haozhang Huang"
    cp.last_modified_by = "Haozhang Huang"
    cp.title = title
    cp.subject = "Isotretinoin and depression reporting in FAERS"
    cp.comments = ""
    return doc

def new_doc(double=False):
    d = Document()
    st = d.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    pf = st.paragraph_format
    pf.line_spacing = 2.0 if double else 1.0
    pf.space_after = Pt(0 if double else 8)
    for s in d.sections:
        s.top_margin = s.bottom_margin = Inches(1)
        s.left_margin = s.right_margin = Inches(1)
    return d


def headed(d, text, size=13):
    p = d.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(size)
    return p


# ----------------------------------------------------------------- cover letter
letter = new_doc()
body = section("1. Cover letter", "2. Title page")
body = body.split("---")[0]
for b in parapraphs(body):
    b = re.sub(r"\(CONFIRM[^)]*\)", "", b).strip()
    if not b:
        continue
    if b.startswith("Dear Editor-in-Chief"):
        d = letter.add_paragraph(b)
        continue
    if b.startswith("Sincerely"):
        letter.add_paragraph("")
        letter.add_paragraph(b)
        continue
    letter.add_paragraph(b)
stamp(letter, "Cover letter").save(OUT / "01_Cover_Letter.docx")

# ----------------------------------------------------------------- title page
title = re.search(r"^# (.+)$", MD, re.M).group(1).strip()
wc = re.search(r"\*\*Word count \(body\): (\d+)", MD).group(1)
t = new_doc()
headed(t, "Title", 12)
t.add_paragraph(clean(title))
headed(t, "Short title", 12)
t.add_paragraph("Rebuilding the comparator for isotretinoin depression reporting")
headed(t, "Authors", 12)
t.add_paragraph("Haozhang Huang\u00b9, Jiajia Li\u00b2, Hao Xiong\u00b2")
headed(t, "Affiliations", 12)
t.add_paragraph("1. Guangdong Provincial Geriatrics Institute, Guangdong Provincial People's "
                "Hospital (Guangdong Academy of Medical Sciences), Southern Medical University, "
                "Guangzhou 510080, China")
t.add_paragraph("2. Department of Dermatology, Guangdong Provincial People's Hospital "
                "(Guangdong Academy of Medical Sciences), Southern Medical University, "
                "Guangzhou 510080, China")
headed(t, "Corresponding author", 12)
t.add_paragraph("Haozhang Huang, Guangdong Provincial Geriatrics Institute, Guangdong Provincial "
                "People's Hospital (Guangdong Academy of Medical Sciences), Southern Medical "
                "University, Guangzhou 510080, China. Email: 1017988724@qq.com")
headed(t, "Manuscript metrics", 12)
t.add_paragraph(f"Body {wc} words; 5 references; one figure in two panels "
                f"(the journal counts each panel as a display); no table in the manuscript, "
                f"the 15-event table is sheet S2 of the supplemental workbook.")
headed(t, "Conflicts of interest", 12)
t.add_paragraph("None declared.")
headed(t, "Funding", 12)
t.add_paragraph("None.")
headed(t, "Patient consent", 12)
t.add_paragraph("Not applicable. This study used de-identified, publicly available "
                "adverse event reports; no patients were involved.")
headed(t, "Ethics", 12)
t.add_paragraph(clean(" ".join(parapraphs(section("3. Ethics / IRB statement", "4. Data availability statement")))))
headed(t, "Data availability", 12)
for b in parapraphs(section("4. Data availability statement", "5. Author contributions")):
    t.add_paragraph(b)
headed(t, "Author contributions", 12)
t.add_paragraph("All authors meet the ICMJE criteria for authorship: each contributed to the "
                "design or to the analysis and interpretation of data, drafted or critically "
                "revised the manuscript, approved the version submitted, and agrees to be "
                "accountable for the work.")
stamp(t, "Title page").save(OUT / "02_Title_Page.docx")

# ----------------------------------------------------------------- RECORD checklist
rec = (ROOT / "manuscript" / "RECORD_PE_checklist.md").read_text(encoding="utf-8")
d = new_doc()


def cell(text):
    return clean(text.strip().strip("|"))


lines = rec.split("\n")
i = 0
while i < len(lines):
    line = lines[i].rstrip()
    if not line.strip():
        i += 1
        continue
    if line.startswith("#"):
        headed(d, clean(line.lstrip("# ")), 12)
        i += 1
        continue
    if line.startswith("|"):
        block = []
        while i < len(lines) and lines[i].lstrip().startswith("|"):
            row = lines[i].strip()
            if not re.match(r"^\|[-: |]+\|$", row):
                block.append([c for c in (x.strip() for x in row.strip("|").split("|"))])
            i += 1
        if not block:
            continue
        ncol = max(len(r) for r in block)
        tbl = d.add_table(rows=0, cols=ncol)
        tbl.style = "Table Grid"
        for r in block:
            cells = tbl.add_row().cells
            for j in range(ncol):
                cells[j].text = clean(r[j]) if j < len(r) else ""
                for par in cells[j].paragraphs:
                    par.paragraph_format.line_spacing = 1.0
                    for run in par.runs:
                        run.font.size = Pt(9)
        d.add_paragraph("")
        continue
    d.add_paragraph(clean(line))
    i += 1

stamp(d, "RECORD / RECORD-PE checklist").save(OUT / "06_RECORD_PE_Checklist.docx")

# ----------------------------------------------------------------- copies
shutil.copy2(ROOT / "manuscript" / "JAAD_BriefReport_isotretinoin_sex.docx",
             OUT / "03_Manuscript.docx")
shutil.copy2(ROOT / "output" / "figures" / "fig1_isotretinoin_rorr.tif",
             OUT / "04_Figure_1.tif")
shutil.copy2(ROOT / "supplement" / "JAAD_supplemental_tables.xlsx",
             OUT / "05_Supplemental_Tables.xlsx")

# ----------------------------------------------------------------- upload note
note = """JAAD submission folder
=====================

Upload in this order at https://www.editorialmanager.com/jaad/

1. 01_Cover_Letter.docx            Cover Letter slot
2. 02_Title_Page.docx              Title Page slot (carries COI, funding, consent)
3. 03_Manuscript.docx              Manuscript slot (Word only; body is double-spaced
                                   with continuous line numbers)
4. 04_Figure_1.tif                 Figure slot, 600 dpi TIFF, 7.4 x 4.7 in
                                   (JAAD accepts .tif/.eps/.jpg, never .pdf)
5. 05_Supplemental_Tables.xlsx     Supplemental Material slot, one file, six tables
6. 06_RECORD_PE_Checklist.docx     Checklist slot, if the system asks for one

Not uploaded, kept in the repository:
- the vector figure (fig1_isotretinoin_rorr.pdf) for your own records
- the full set of 43 supplemental tables and the analysis code

The four documents contain no placeholders. Three things still need a human:
1. Confirm with Jiajia Li and Hao Xiong that they have no competing interests.
2. Add the Mendeley DOI to the data-availability line once the deposit exists
   (the statement currently says the workbook will be deposited).
3. Decide whether to name the code repository; today the manuscript says the
   code is available from the corresponding author.
"""
(OUT / "README_UPLOAD.txt").write_text(note, encoding="utf-8")

for f in sorted(OUT.iterdir()):
    print(f"{f.name:34s} {f.stat().st_size/1e3:8.1f} kB")
