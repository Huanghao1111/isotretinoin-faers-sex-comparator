#!/usr/bin/env python
# ==============================================================================
# Freeze what is being submitted.
#
# The folder was rebuilt six times today and every rebuild changed the hashes.
# A journal that later asks "which version did you send" needs an answer that is
# not "the folder as it stands now". This writes a dated archive and a manifest
# recording the hashes, the word count and the state of the verification suite.
# ==============================================================================
import hashlib
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path

ROOT = Path("C:/Users/lenovo/WorkBuddy/2026-10-17-00-49-45")
SUB = ROOT / "submission"
REL = ROOT / "release"
REL.mkdir(exist_ok=True)

stamp = datetime.now().strftime("%Y-%m-%d")
out = REL / f"JAAD_submission_{stamp}.zip"

suites = {
    "08 numbers in the text": ["Rscript", "code/08_verify_numbers.R"],
    "28 supplement and folder": ["python", "code/28_check_supplement_links.py"],
    "34 retired phrases": ["python", "code/34_check_retired_phrases.py"],
    "35 numbers match their table": ["python", "code/35_check_checklist_numbers.py"],
    "99 figure": ["Rscript", "code/99_check_pdf.R"],
    "23c Word file": ["python", "code/23c_check_docx.py"],
}

lines = [f"JAAD submission, frozen {datetime.now():%Y-%m-%d %H:%M}", ""]
lines.append(f"{'file':34s} {'bytes':>9s}  sha256")
for f in sorted(SUB.iterdir()):
    if f.is_file():
        h = hashlib.sha256(f.read_bytes()).hexdigest()
        lines.append(f"{f.name:34s} {f.stat().st_size:9d}  {h[:32]}")

lines += ["", "verification at freeze:"]
for name, cmd in suites.items():
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    tail = [l for l in r.stdout.splitlines() if "passed" in l or "failures:" in l]
    lines.append(f"  {name:32s} {' '.join(tail[-1:]).strip() or '(no summary)'}")

manifest = "\n".join(lines) + "\n"
(REL / f"MANIFEST_{stamp}.txt").write_text(manifest, encoding="utf-8")

with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(SUB.iterdir()):
        if f.is_file():
            z.write(f, f.name)
    z.writestr("MANIFEST.txt", manifest)

print(manifest)
print(f"[written] {out} ({out.stat().st_size/1e6:.2f} MB)")
