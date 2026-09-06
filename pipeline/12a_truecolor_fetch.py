"""True-colour check windows for a tile: NWIN square windows of SIDE m placed where sites are densest (non-overlapping),
each fetched on the wettest clear scene of its two wettest seasons. Needs only earthengine-api + numpy (gee env).
Writes out/12_truecolor/raw/{TILE}_w{k}_{year}.npy (uint8 RGB) + out/12_truecolor/raw/{TILE}_windows.json.
Run:  TILE=19_39 nice -n 15 $(cat .python_env_gee) scripts/12a_truecolor_fetch.py     (NWIN=6 SIDE=12000)"""
import os, io, csv, json, glob, urllib.request, numpy as np, ee
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); RAW = os.path.join(OUT, "12_truecolor", "raw")
TILE = os.environ.get("TILE", "19_39"); NWIN = int(os.environ.get("NWIN", "6")); SIDE = float(os.environ.get("SIDE", "12000")); SCALE = 10
ee.Initialize(project=open(os.path.join(ROOT, ".gee_project")).read().strip())
rr, cc = (int(v) for v in TILE.split("_")); X0 = (cc - 1) * 100_000 - 4_000_000; Y0 = (rr - 1) * 100_000 - 4_000_000
reg = list(csv.DictReader(open(os.path.join(OUT, f"09_registry_{TILE}.csv"))))
gj = json.load(open(os.path.join(OUT, f"09_sites_{TILE}.geojson"))); geom_of = {f["properties"]["site_id"]: f["geometry"] for f in gj["features"]}; sy = {r["site_id"]: r for r in csv.DictReader(open(os.path.join(OUT, f"09_site_years_{TILE}.csv")))}
years = sorted(int(k[1:]) for k in next(iter(sy.values())) if k.startswith("a"))
# windows: grid cells of SIDE, ranked by site count, greedy non-overlapping
cells = {}
for r in reg:
    k = (int((float(r["x3413"]) - X0) // SIDE), int((float(r["y3413"]) - Y0) // SIDE)); cells.setdefault(k, []).append(r)
chosen = []
for k, rs in sorted(cells.items(), key=lambda kv: -sum(float(r['max_year_area_km2']) for r in kv[1] if r.get('ice_marginal', 'False') != 'True')):  # rank windows by peak lake area, not site count
    if len(chosen) >= NWIN: break
    chosen.append((k, rs))
def prep(img):
    b = img.select(["B2", "B3", "B4", "B11"]).divide(1e4); ndsi = b.normalizedDifference(["B3", "B11"]); margin = ndsi.lt(0.85).And(b.select("B2").lt(0.4))
    valid = b.select("B2").mask().And(b.select("B4").mask()).And(margin.Not()); ndwi = b.normalizedDifference(["B2", "B4"])
    return ee.Image.cat(valid.rename("obs"), ndwi.gt(0.5).And(valid).rename("w")).unmask(0)
DROP = {y: json.load(open(f)).get("dropped_toa", []) for y in years for f in glob.glob(os.path.join(OUT, f"08_gee_task_{TILE}_{y}.json"))}
def best_scene(region, year, sites_geom):
    col = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").filterDate(f"{year}-07-01", f"{year}-09-16").filterBounds(region).filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 60))
    if DROP.get(year): col = col.filter(ee.Filter.inList("system:index", DROP[year]).Not())
    cp = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY").filterDate(f"{year}-07-01", f"{year}-09-16").filterBounds(region)
    col = ee.ImageCollection(ee.Join.saveFirst("cloud").apply(col, cp, ee.Filter.equals(leftField="system:index", rightField="system:index")))
    def stat(img):
        pr = prep(img); s = pr.select("obs").reduceRegion(ee.Reducer.sum(), region, scale=60, crs="EPSG:3413", maxPixels=1e7)
        c = ee.Image(img.get("cloud")).select("probability").reduceRegion(ee.Reducer.mean(), region, scale=120, crs="EPSG:3413", maxPixels=1e7)
        w = pr.select("w").reduceRegion(ee.Reducer.sum(), sites_geom, scale=30, crs="EPSG:3413", maxPixels=1e7); return img.set({"obs": s.get("obs"), "w": w.get("w"), "cloud": c.get("probability")})
    col = col.map(stat); mx = ee.Number(col.aggregate_max("obs")); clear = col.filter(ee.Filter.gte("obs", mx.multiply(0.97))).filter(ee.Filter.lt("cloud", 20))
    return ee.Image(ee.Algorithms.If(clear.size().gt(0), clear.sort("w", False).first(), col.filter(ee.Filter.gte("obs", mx.multiply(0.97))).sort("cloud").first()))
def rgb(img, wx0, wy1, wpx):
    im = img.select(["B4", "B3", "B2"]).divide(1e4).divide(0.7).clamp(0, 1).multiply(255).toUint8()
    url = im.getDownloadURL({"crs": "EPSG:3413", "crs_transform": [SCALE, 0, wx0, 0, -SCALE, wy1], "dimensions": f"{wpx}x{wpx}", "format": "NPY", "bands": ["B4", "B3", "B2"]})
    a = np.load(io.BytesIO(urllib.request.urlopen(url, timeout=600).read())); return np.dstack([a["B4"], a["B3"], a["B2"]])
wins = []; wpx = int(SIDE / SCALE)
for k, (cell, rs) in enumerate(chosen):
    wx0 = X0 + cell[0] * SIDE; wy0 = Y0 + cell[1] * SIDE; wx1, wy1 = wx0 + SIDE, wy0 + SIDE
    region = ee.Geometry.Rectangle([wx0, wy0, wx1, wy1], proj="EPSG:3413", evenOdd=False)
    keep = [r for r in rs if r.get("ice_marginal", "False") != "True"] or rs
    sites_geom = ee.FeatureCollection([ee.Feature(ee.Geometry(geom_of[r["site_id"]])) for r in keep]).geometry()
    area = {y: sum(float(sy[r["site_id"]][f"a{y}"]) for r in keep) for y in years}; seas = sorted(sorted(years, key=lambda y: -area[y])[:2]); dates = {}
    for y in seas:
        f = os.path.join(RAW, f"{TILE}_w{k}_{y}.npy")
        if os.path.exists(f): dates[y] = "cached"; continue
        try:
            img = best_scene(region, y, sites_geom); dates[y] = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd").getInfo(); np.save(f, rgb(img, wx0, wy1, wpx))
        except Exception as e: dates[y] = f"ERROR {str(e)[:80]}"
    wins.append(dict(k=k, window=[wx0, wy0, wx1, wy1], n_sites=len(rs), seasons=seas, dates=dates)); print(k, len(rs), "sites", dates, flush=True)
    json.dump(wins, open(os.path.join(RAW, f"{TILE}_windows.json"), "w"), indent=1)
print("done")
