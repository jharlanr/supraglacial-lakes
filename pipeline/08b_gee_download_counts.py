"""Download the s2counts asset made by 08a as numpy arrays, in windows small enough for
getDownloadURL (2500 x 2500 px x 4 uint8 bands = 25 MB each; 16 windows for a 100 km tile at 10 m).
Writes out/08_s2counts_{TILE}_{YEAR}_{band}.npy (uint8, 10000 x 10000, row 0 = north) + _meta.json.

Run:  $(cat .python_env_gee) scripts/08b_gee_download_counts.py        (TILE=19_39 YEAR=2019)
"""
import os, io, json, time, urllib.request
import numpy as np, ee
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out")
TILE = os.environ.get("TILE", "19_39"); YEAR = int(os.environ.get("YEAR", "2019")); WIN = int(os.environ.get("WIN", "2500"))
info = json.load(open(os.path.join(OUT, f"08_gee_task_{TILE}_{YEAR}.json")))
ee.Initialize(project=open(os.path.join(ROOT, ".gee_project")).read().strip())
st = ee.data.getTaskStatus(info["task_id"])[0]; assert st["state"] == "COMPLETED", st
img = ee.Image(info["asset"]); bands = img.bandNames().getInfo(); print("bands:", bands)
x0, y0, x1, y1 = info["x0"], info["y0"], info["x1"], info["y1"]
N = (x1 - x0) // 10; arrs = {b: np.zeros((N, N), np.uint8) for b in bands}
for r0 in range(0, N, WIN):
    for c0 in range(0, N, WIN):
        wx0, wy1 = x0 + c0 * 10, y1 - r0 * 10; h = min(WIN, N - r0); w = min(WIN, N - c0)
        for attempt in range(4):
            try:
                url = img.getDownloadURL({"crs": "EPSG:3413", "crs_transform": [10, 0, wx0, 0, -10, wy1],
                                          "dimensions": f"{w}x{h}", "format": "NPY", "bands": bands})
                raw = urllib.request.urlopen(url, timeout=600).read(); a = np.load(io.BytesIO(raw))
                for b in bands: arrs[b][r0:r0 + h, c0:c0 + w] = a[b]
                print(f"window r{r0} c{c0}: {len(raw)/1e6:.1f} MB"); break
            except Exception as e:
                print("retry", attempt, type(e).__name__, str(e)[:120]); time.sleep(10 * (attempt + 1))
        else: raise SystemExit("window failed")
for b in bands: np.save(os.path.join(OUT, f"08_s2counts_{TILE}_{YEAR}_{b}.npy"), arrs[b])
meta = dict(info, bands=bands, transform=[10, 0, x0, 0, -10, y1], shape=[N, N], crs="EPSG:3413",
            downloaded=time.strftime("%Y-%m-%d %H:%M %Z"), n_obs_median=int(np.median(arrs["n_obs"])))
json.dump(meta, open(os.path.join(OUT, f"08_s2counts_{TILE}_{YEAR}_meta.json"), "w"), indent=1)
print("saved", {b: f"{arrs[b].max()} max" for b in bands}, "| n_obs median", meta["n_obs_median"])
