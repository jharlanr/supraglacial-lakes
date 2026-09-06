"""Fetch Sentinel-2 true-colour windows for the closing-accuracy cases (out/11_cases_{TILE}.csv from 11a): for each case's
NSEAS wettest seasons, the scene with the most water pixels among scenes covering >= 95 % of the window.
Writes out/11_chips/raw/{case_id}_{year}.npy (uint8 RGB) and out/11_chips/raw/{case_id}.json (window, seasons, dates).
Needs only earthengine-api + numpy.   Run:  TILE=19_39 nice -n 15 $(cat .python_env_gee) scripts/11b_closing_fetch.py   (CASES=n limits)"""
import os, io, csv, json, glob, urllib.request, numpy as np, ee
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); RAW = os.path.join(OUT, "11_chips", "raw")
TILE = os.environ.get("TILE", "19_39"); LIMIT = int(os.environ.get("CASES", "0")); NSEAS = 4
ee.Initialize(project=open(os.path.join(ROOT, ".gee_project")).read().strip())
cases = list(csv.DictReader(open(os.path.join(OUT, f"11_cases_{TILE}.csv")))); cases = cases[:LIMIT] if LIMIT else cases
sy = {r["site_id"]: r for r in csv.DictReader(open(os.path.join(OUT, f"09_site_years_{TILE}.csv")))}
years = sorted(int(os.path.basename(f)[len(f"08_s2counts_{TILE}_"):].split("_")[0]) for f in glob.glob(os.path.join(OUT, f"08_s2counts_{TILE}_*_meta.json")))
def prep(img):
    b = img.select(["B2", "B3", "B4", "B11"]).divide(1e4); ndsi = b.normalizedDifference(["B3", "B11"]); margin = ndsi.lt(0.85).And(b.select("B2").lt(0.4))
    valid = b.select("B2").mask().And(b.select("B4").mask()).And(margin.Not()); ndwi = b.normalizedDifference(["B2", "B4"])
    return ee.Image.cat(valid.rename("obs"), ndwi.gt(0.5).And(valid).rename("w")).unmask(0)
DROP = {y: json.load(open(f)).get("dropped_toa", []) for y in years for f in glob.glob(os.path.join(OUT, f"08_gee_task_{TILE}_{y}.json"))}
def best_scene(region, year):
    col = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").filterDate(f"{year}-06-01", f"{year}-10-01").filterBounds(region)
    if DROP.get(year): col = col.filter(ee.Filter.inList("system:index", DROP[year]).Not())
    def stat(img):
        s = prep(img).reduceRegion(ee.Reducer.sum(), region, scale=30, crs="EPSG:3413", maxPixels=1e7); return img.set({"obs": s.get("obs"), "w": s.get("w")})
    col = col.map(stat); mx = ee.Number(col.aggregate_max("obs"))
    return ee.Image(col.filter(ee.Filter.gte("obs", mx.multiply(0.95))).sort("w", False).first())
def rgb(img, wx0, wy1, wpx, scale):
    im = img.select(["B4", "B3", "B2"]).divide(1e4).divide(0.7).clamp(0, 1).multiply(255).toUint8()
    url = im.getDownloadURL({"crs": "EPSG:3413", "crs_transform": [scale, 0, wx0, 0, -scale, wy1], "dimensions": f"{wpx}x{wpx}", "format": "NPY", "bands": ["B4", "B3", "B2"]})
    a = np.load(io.BytesIO(urllib.request.urlopen(url, timeout=300).read())); return np.dstack([a["B4"], a["B3"], a["B2"]])
for c in cases:
    jp = os.path.join(RAW, f"{c['case_id']}.json")
    if os.path.exists(jp): continue
    xmin, xmax, ymin, ymax = (float(c[k]) for k in ("xmin", "xmax", "ymin", "ymax"))
    side = max(1200.0, max(xmax - xmin, ymax - ymin) + 600); cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    scale = 10 if side <= 4000 else 20; wpx = int(round(side / scale)); wx0 = round((cx - side / 2) / 10) * 10; wy1 = round((cy + side / 2) / 10) * 10; wx1, wy0 = wx0 + wpx * scale, wy1 - wpx * scale
    region = ee.Geometry.Rectangle([wx0, wy0, wx1, wy1], proj="EPSG:3413", evenOdd=False)
    sids = c["site_ids"].split(";"); area = np.sum([[float(sy[s][f"a{y}"]) for y in years] for s in sids], axis=0)
    seas = sorted(years[i] for i in np.argsort(-area, kind="stable")[:NSEAS]); dates = {}
    for y in seas:
        try:
            img = best_scene(region, y); dates[y] = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd").getInfo()
            np.save(os.path.join(RAW, f"{c['case_id']}_{y}.npy"), rgb(img, wx0, wy1, wpx, scale))
        except Exception as e: dates[y] = f"ERROR {str(e)[:80]}"
    json.dump(dict(case_id=c["case_id"], pool=c["pool"], site_ids=sids, n_pieces=int(float(c["n_pieces"])), gap_m=None if c["gap_m"] in ("", "nan") else int(float(c["gap_m"])),
                   window=[wx0, wy0, wx1, wy1], scale=scale, wpx=wpx, seasons=seas, dates=dates), open(jp, "w"))
    print(c["case_id"], dates, f"{side:.0f} m", flush=True)
print("fetch done")
