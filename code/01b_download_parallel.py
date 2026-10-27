"""
Concurrent FAERS quarterly downloader.
Writes to *.part then renames, so a partial file is never mistaken for complete.
"""
import os
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

RAW = r"C:\Users\lenovo\WorkBuddy\2026-10-17-00-49-45\data\raw"
os.makedirs(RAW, exist_ok=True)

BASE = "https://fis.fda.gov/content/Exports/faers_ascii_{y}q{q}.zip"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
MIN_SIZE = 500_000
WORKERS = 8


def quarters(start=(2012, 4), end=(2026, 4)):
    y, q = start
    while (y, q) <= end:
        yield y, q
        q += 1
        if q > 4:
            q = 1
            y += 1


def fetch(y, q, tries=3):
    tag = f"{y}q{q}"
    dest = os.path.join(RAW, f"faers_ascii_{tag}.zip")
    part = dest + ".part"
    if os.path.exists(dest) and os.path.getsize(dest) >= MIN_SIZE:
        return tag, "skip", 0
    for attempt in range(tries):
        try:
            req = urllib.request.Request(BASE.format(y=y, q=q),
                                         headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            if len(data) < MIN_SIZE:
                return tag, f"notfound({len(data)}B)", 0
            with open(part, "wb") as f:
                f.write(data)
            os.replace(part, dest)
            return tag, "ok", len(data)
        except Exception as e:
            if attempt == tries - 1:
                return tag, f"fail:{e}", 0
            time.sleep(3)
    return tag, "fail", 0


def main():
    tasks = list(quarters())
    ok = skip = miss = 0
    failed = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(fetch, y, q) for y, q in tasks]
        for fu in as_completed(futs):
            tag, status, size = fu.result()
            if status == "ok":
                ok += 1
                print(f"[OK]   {tag}  {size/1e6:.1f} MB", flush=True)
            elif status == "skip":
                skip += 1
            elif status.startswith("notfound"):
                miss += 1
            else:
                failed.append(tag)
                print(f"[ERR]  {tag}  {status}", flush=True)
    print("\n==== SUMMARY ====")
    print(f"downloaded={ok} skipped={skip} unavailable={miss} failed={len(failed)}")
    if failed:
        print("failed:", " ".join(failed))


if __name__ == "__main__":
    main()
