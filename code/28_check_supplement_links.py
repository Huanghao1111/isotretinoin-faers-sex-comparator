#!/usr/bin/env python
# ==============================================================================
# Checks the supplemental numbering: every deposited table is described, every
# reference points at a table that exists, and the numbering has no holes.
#
# This is here because the package has drifted twice during writing - a table was
# described as "rejected" after it had been superseded, and the S-range quoted in
# three documents lagged the files by several entries. Both are invisible to
# proof-reading and obvious to a script.
#
# Exit status is non-zero on any failure.
# ==============================================================================
import re
import sys
from pathlib import Path

ROOT = Path("C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45")
SUP = ROOT / "supplement"
DOCS = [ROOT / "supplement" / "README.md"] + sorted((ROOT / "manuscript").glob("*.md"))
# The three review reports describe the package as it stood when they were run;
# their quoted S-range is part of that record and is not updated retrospectively.
LIVE = [d for d in DOCS if not d.name.startswith("99")]

tables = {}
dupes = []
for f in sorted(SUP.glob("TableS*_*")):
    m = re.match(r"TableS(\d+)_", f.name)
    if not m:
        continue
    n = int(m.group(1))
    if n in tables:
        dupes.append((n, tables[n].name, f.name))
    tables[n] = f

ok = fail = 0


def check(label, cond, detail=""):
    global ok, fail
    if cond:
        ok += 1
        print(f"[PASS] {label:<52} {detail}")
    else:
        fail += 1
        print(f"[FAIL] {label:<52} {detail}")


nums = sorted(tables)
check("tables deposited", len(nums) > 0, f"{len(nums)} files, S{nums[0]}-S{nums[-1]}")
check("numbering contiguous", nums == list(range(1, max(nums) + 1)),
      "no gaps" if nums == list(range(1, max(nums) + 1)) else f"gaps: {sorted(set(range(1, max(nums)+1)) - set(nums))}")
check("no duplicate numbers", not dupes, str(dupes) if dupes else "unique")

# every reference in every document must point at a deposited table
refs = {}
for d in DOCS:
    text = d.read_text(encoding="utf-8")
    for m in re.finditer(r"TableS(\d+)", text):
        refs.setdefault(int(m.group(1)), set()).add(d.name)

missing = sorted(n for n in refs if n not in tables)
check("all references resolve", not missing,
      "every TableSxx has a file" if not missing else f"dangling: {missing}")

# every deposited table must be described in the README
readme = (SUP / "README.md").read_text(encoding="utf-8")
undescribed = [n for n in nums if f"TableS{n}_" not in readme and f"TableS{n} " not in readme
               and f"TableS{n}." not in readme and f"TableS{n}\n" not in readme]
check("every table described in README", not undescribed,
      "all described" if not undescribed else f"undescribed: {undescribed}")

# the S-range quoted in prose must match the files
patterns = [
    (r"TableS1[-\u2013]S(\d+)", "TableS1-Sxx"),
    (r"S1 to S(\d+)", "S1 to Sxx"),
]
for d in LIVE:
    text = d.read_text(encoding="utf-8")
    for pat, label in patterns:
        for m in re.finditer(pat, text):
            claimed = int(m.group(1))
            check(f"{d.name}: {label} range", claimed == max(nums),
                  f"claimed S{claimed}, files end at S{max(nums)}")

# and the tables referenced from the manuscript itself
# the Word file and the repository README are uploaded/deposited too
docx = ROOT / "manuscript" / "JAAD_BriefReport_isotretinoin_sex.docx"
if docx.exists():
    import zipfile, re as _re
    with zipfile.ZipFile(docx) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    text = _re.sub(r"<[^>]+>", "", xml)
    for m in re.finditer(r"TableS1[-\u2013]S(\d+)", text):
        check("submission docx: TableS1-Sxx range", int(m.group(1)) == max(nums),
              f"claimed S{m.group(1)}, files end at S{max(nums)}")
root_readme = ROOT / "README.md"
if root_readme.exists():
    text = root_readme.read_text(encoding="utf-8")
    for m in re.finditer(r"Tables? S1[-\u2013]S(\d+)", text):
        check("repository README: S-range", int(m.group(1)) == max(nums),
              f"claimed S{m.group(1)}, files end at S{max(nums)}")

# the single file that is actually deposited
wb = SUP / "JAAD_supplemental_tables.xlsx"
if wb.exists():
    import zipfile, re as _re
    with zipfile.ZipFile(wb) as z:
        wbxml = z.read("xl/workbook.xml").decode("utf-8")
    sheet_names = _re.findall(r'name="(S\d+)"', wbxml)
    expected = [f"S{i}" for i in range(1, 7)]
    check("deposit workbook: six tables", sheet_names == expected, f"sheets {sheet_names}")
    src = (SUP / "JAAD_supplemental_tables.xlsx").stat().st_size / 1e6
    check("deposit workbook: file size", src < 5, f"{src:.2f} MB")

ms = (ROOT / "manuscript" / "JAAD_BriefReport_isotretinoin_sex.md").read_text(encoding="utf-8")
ms_refs = sorted({int(x) for x in re.findall(r"TableS(\d+)", ms)})
check("manuscript references deposited tables",
      all(n in tables for n in ms_refs), f"refs: {ms_refs}" if ms_refs else "none cited")

# the upload folder must be complete and current
sub = ROOT / "submission"
if sub.exists():
    need = ["01_Cover_Letter.docx", "02_Title_Page.docx", "03_Manuscript.docx",
            "04_Figure_1.tif", "05_Supplemental_Tables.xlsx",
            "06_RECORD_PE_Checklist.docx", "README_UPLOAD.txt"]
    missing = [n for n in need if not (sub / n).exists()]
    check("submission folder complete", not missing, str(missing) if missing else "7 files")
    same = (sub / "03_Manuscript.docx").read_bytes() == \
           (ROOT / "manuscript" / "JAAD_BriefReport_isotretinoin_sex.docx").read_bytes()
    check("submission manuscript is the current build", same,
          "identical" if same else "STALE - rerun 31")
    same = (sub / "05_Supplemental_Tables.xlsx").read_bytes() == \
           (SUP / "JAAD_supplemental_tables.xlsx").read_bytes()
    check("submission supplement is the current workbook", same,
          "identical" if same else "STALE - rerun 31")

print(f"\n===== {ok} passed, {fail} failed =====")
sys.exit(1 if fail else 0)
