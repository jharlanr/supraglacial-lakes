"""Submit one Earth Engine batch task: per-pixel Sentinel-2 water-detection counts for one
ArcticDEM tile (EPSG:3413, 100 km) and one melt season, at 10 m, using EVERY scene (no
scene-level cloud filter). Bands (uint8):
  n_obs        scenes with valid (unmasked, non-margin) data at the pixel      [TOA collection]
  n_w50_toa    scenes with NDWI(B2,B4) > 0.5                                    [TOA]
  n_w30_toa    scenes with NDWI > 0.3                                           [TOA]
  n_w50_sr     scenes with NDWI > 0.5 on the surface-reflectance collection     [SR, where it exists]
Margin/rock/seawater mask per scene: NDSI = (B3-B11)/(B3+B11) < 0.85 AND B2 < 0.4 (Moussavi 2020,
as in Dunmire 2021). The asset lands in the project; 08b downloads it.

Run:  $(cat .python_env_gee) scripts/08a_gee_s2_counts.py        (TILE=19_39 YEAR=2019 by default)
"""
import os, ee, json, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT = open(os.path.join(ROOT, ".gee_project")).read().strip()
TILE = os.environ.get("TILE", "19_39"); YEAR = int(os.environ.get("YEAR", "2019"))
T0, T1 = os.environ.get("T0", f"{YEAR}-06-01"), os.environ.get("T1", f"{YEAR}-10-01")
ee.Initialize(project=PROJECT)

rr, cc = (int(v) for v in TILE.split("_"))
x0 = -4_000_000 + (cc - 1) * 100_000; y0 = -4_000_000 + (rr - 1) * 100_000
x1, y1 = x0 + 100_000, y0 + 100_000
region = ee.Geometry.Rectangle([x0, y0, x1, y1], proj="EPSG:3413", evenOdd=False)

def prep(img, scale=1e4):
    b = img.select(["B2", "B3", "B4", "B11"]).divide(scale)
    ndsi = b.normalizedDifference(["B3", "B11"])
    margin = ndsi.lt(0.85).And(b.select("B2").lt(0.4))          # rock / seawater -> excluded
    valid = b.select("B2").mask().And(b.select("B4").mask()).focalMin(EDGE_PX, "square", "pixels").And(margin.Not())  # drop EDGE_PX swath-edge pixels (line artefacts)
    ndwi = b.normalizedDifference(["B2", "B4"])
    return ee.Image.cat(valid.rename("obs"), ndwi.gt(0.5).And(valid).rename("w50"), ndwi.gt(0.3).And(valid).rename("w30")).unmask(0)

# Scene QA (added 2026-09-05 after the 2024-08-31 T22WEC striped scene put 89 km2 of diagonal bands into a season):
# water area per scene over the tile at 60 m; scenes above QA_FACTOR x the season's 95th percentile are dropped.
QA_FACTOR = float(os.environ.get("QA_FACTOR", "3")); EDGE_PX = int(os.environ.get("EDGE_PX", "3"))
def with_wpx(img):
    b = img.select(["B2", "B4"]).divide(1e4); w = b.normalizedDifference(["B2", "B4"]).gt(0.5).rename("w")
    return img.set("wpx", w.reduceRegion(ee.Reducer.sum(), region, scale=60, maxPixels=1e9).get("w"))
def qa(col):
    col = col.map(with_wpx); vals = col.aggregate_array("wpx")
    # anchor on the top of the distribution: the median of the ten largest scene water areas. A season's 95th
    # percentile is near zero in low-melt years (most scenes hold no water) and would reject the peak scenes.
    top10 = vals.sort().reverse().slice(0, 10); anchor = ee.Number(top10.reduce(ee.Reducer.median()))
    cut = anchor.multiply(QA_FACTOR).max(1)
    kept = col.filter(ee.Filter.lte("wpx", cut)); dropped = col.filter(ee.Filter.gt("wpx", cut))
    return kept, dropped, cut
toa_all = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").filterDate(T0, T1).filterBounds(region)
sr_all = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterDate(T0, T1).filterBounds(region)
toa_k, toa_d, toa_cut = qa(toa_all); sr_k, sr_d, sr_cut = qa(sr_all)
dropped_toa = toa_d.aggregate_array("system:index").getInfo(); dropped_sr = sr_d.aggregate_array("system:index").getInfo()
toa = toa_k.map(prep); sr = sr_k.map(prep)
n_toa, n_sr = toa.size().getInfo(), sr.size().getInfo()
print(f"tile {TILE} {YEAR}: {n_toa} TOA scenes kept, {len(dropped_toa)} dropped by QA (cut {toa_cut.getInfo():.0f} px at 60 m): {dropped_toa}; {n_sr} SR kept, {len(dropped_sr)} dropped")

counts = ee.Image.cat(
    toa.select("obs").sum().rename("n_obs"),
    toa.select("w50").sum().rename("n_w50_toa"),
    toa.select("w30").sum().rename("n_w30_toa"),
    sr.select("w50").sum().rename("n_w50_sr"),
).toUint8().set({"tile": TILE, "year": YEAR, "t0": T0, "t1": T1, "n_toa_scenes": n_toa, "n_sr_scenes": n_sr,
                 "recipe": "every scene passing scene-QA; margin NDSI<0.85&B2<0.4; NDWI(B2,B4)>0.5|0.3; no cloud filter", "qa_factor": QA_FACTOR, "dropped_toa": ",".join(dropped_toa), "dropped_sr": ",".join(dropped_sr)})

asset = f"projects/{PROJECT}/assets/s2counts_{TILE}_{YEAR}"
if os.environ.get("OVERWRITE") == "1":
    try: ee.data.deleteAsset(asset); print("deleted existing asset")
    except Exception as e: print("no existing asset to delete:", str(e)[:80])
task = ee.batch.Export.image.toAsset(image=counts, description=f"s2counts_{TILE}_{YEAR}", assetId=asset,
                                     region=region, crs="EPSG:3413", crsTransform=[10, 0, x0, 0, -10, y1],
                                     maxPixels=2e9, pyramidingPolicy={".default": "sample"})
task.start()
info = {"task_id": task.id, "asset": asset, "tile": TILE, "year": YEAR, "x0": x0, "y0": y0, "x1": x1, "y1": y1,
        "n_toa_scenes": n_toa, "n_sr_scenes": n_sr, "dropped_toa": dropped_toa, "dropped_sr": dropped_sr, "submitted": time.strftime("%Y-%m-%d %H:%M %Z")}
json.dump(info, open(os.path.join(ROOT, "out", f"08_gee_task_{TILE}_{YEAR}.json"), "w"), indent=1)
print("submitted task", task.id, "->", asset)
