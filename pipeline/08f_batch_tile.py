"""Submit the ten-season Earth Engine exports for one tile (08a), then poll and download each on completion (08b).
Run:  TILE=29_45 $(cat .python_env_gee) scripts/08f_batch_tile.py      (background; log to out/08_batch_{TILE}.log)
Skips seasons whose out/08_s2counts_{TILE}_{YEAR}_meta.json already exists; resubmits none (delete the task json to force)."""
import os, sys, json, time, subprocess, datetime
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "29_45"); YEARS = [int(y) for y in os.environ.get("YEARS", "2016 2017 2018 2019 2020 2021 2022 2023 2024 2025").split()]
PY = sys.executable
def log(msg): print(datetime.datetime.now().strftime("%H:%M"), msg, flush=True)
pending = {}
for y in YEARS:
    if os.path.exists(os.path.join(OUT, f"08_s2counts_{TILE}_{y}_meta.json")): log(f"{y}: already downloaded"); continue
    tj = os.path.join(OUT, f"08_gee_task_{TILE}_{y}.json")
    if not os.path.exists(tj):
        r = subprocess.run([PY, os.path.join(ROOT, "scripts", "08a_gee_s2_counts.py")], env={**os.environ, "TILE": TILE, "YEAR": str(y)}, capture_output=True, text=True)
        log(f"{y}: submit rc={r.returncode} {r.stdout.strip()[-300:]} {r.stderr.strip()[-300:]}")
    pending[y] = json.load(open(tj))["task_id"]
import ee; ee.Initialize(project=open(os.path.join(ROOT, ".gee_project")).read().strip())
while pending:
    time.sleep(120)
    for y, tid in list(pending.items()):
        st = ee.data.getTaskStatus(tid)[0]["state"]
        if st == "COMPLETED":
            r = subprocess.run([PY, os.path.join(ROOT, "scripts", "08b_gee_download_counts.py")], env={**os.environ, "TILE": TILE, "YEAR": str(y)}, capture_output=True, text=True)
            log(f"{y}: COMPLETED, download rc={r.returncode} {r.stdout.strip()[-200:]}"); pending.pop(y)
        elif st in ("FAILED", "CANCELLED"): log(f"{y}: {st} {ee.data.getTaskStatus(tid)[0].get('error_message','')}"); pending.pop(y)
    log("waiting on " + " ".join(str(y) for y in pending))
log("batch done")
