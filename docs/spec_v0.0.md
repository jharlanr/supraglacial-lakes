# Greenland lake sites — specification v0.0 (draft, 2026-09-05 Pacific)

*Status: design agreed in conversation with Josh on 2026-09-05; nothing built yet. Every number marked
"default" is a knob to be tested before the registry is frozen. Precedents are in
`notes/lit_lake_definition_and_identity.md`.*

## 1. Definitions

**Outline (per year).** The seasonal maximum extent of a lake in one melt season. Not a reproduction of
Dunmire 2021 (Josh, 2026-09-05): native 10 m, no scene-level cloud filter.
- Sentinel-2 L2A scenes over the ice sheet, melt season (default June–September), **every scene**, no
  scene-level cloud filter: clouds are bright and never pass the water threshold, and a scene filter discards
  clear pixels because of cloud elsewhere in the granule. Ice-margin / rock / seawater mask kept (NDSI < 0.85
  and blue < 0.4, Moussavi 2020) — it is not a cloud mask.
- Grid 10 m, native. Landsat back-fill outlines are 30 m; membership by overlap does not care.
- Water pixel: NDWI = (blue − red)/(blue + red) > 0.5 (Dunmire 2021, from Miles 2017). High on purpose:
  in a union, a false positive is permanent and a false negative is filled by another date.
- Scene QA (added 2026-09-05 after one striped acquisition, 2024-08-31, put 89 km² of diagonal bands into a season):
  per scene, the water area over the tile at 60 m; a scene above 3 × the median of the season's ten largest scene
  water areas is dropped. Anchoring on the top of the distribution matters: a season's 95th percentile is near zero
  in low-melt years and would reject the peak scenes (2018: peak 7 × p95; 2025: 11 ×). Over ten seasons of tile
  19_39 the rule removes only the four granules of that one acquisition (`out/09_scene_qa_scan_19_39.json`).
- Persistence: **k = 1** — a pixel is water for the year if it passes the threshold in any scene (Dunmire's own
  rule). Tested 2026-09-05 on tile 19_39, 2019 (`out/08_outline_knobs_19_39_2019.txt`): k = 2 loses 9 of 217
  Dunmire lakes, k = 3 loses 28, because a lake's outline is water in few scenes (median 11 detections in ~97
  observations). Cloud shadow at 0.5 turned out not to produce ≥ 0.05 km² components (no extra component was a
  single-shot detection), so no scene cloud filter and no persistence count are needed.
- Per-year mask = pixels that met the rule → closing to fill floating ice and bridge lid gaps (default disk radius
  5 px = 50 m at 10 m; 0 / 5 / 10 px tested, 5 best) → **no opening** (no river chaining at 0.5: at most 2 Dunmire
  lakes under one component in every configuration; the opening cost 3–9 lakes) → connected components → keep
  area ≥ 0.05 km² (= 500 px) → **fill filter**: keep only components in which ≥ 50 % of the pixels were detected
  water (lakes: 94 % median; slush and crevasse fields bridged by the closing: 49 %). The fill filter halves the
  non-lake components at no cost in Dunmire lakes.
- Union of seasons (tested with 2018+2019, `out/08_multiyear_sites_19_39.txt`): 2019 recall rises from 191 to 203 of 217
  Dunmire lakes and 138 of 139 for 2018 with the per-year recipe + a 150 m site closing; 12 of the 26 single-season misses
  come back; no chaining (≤ 2 lakes per site); the crescent becomes one site from the union alone. A 'seen in ≥ 2 years'
  site rule cannot be judged with two seasons.
- DEM coverage (Friday's overlay, four CW tiles, 1 538 lakes): 9 % of 2018 and 12 % of 2019 Dunmire lakes lie entirely
  outside any 1 m-deep 32 m depression; 27 % of lakes < 0.1 km², 16 % of 0.1–0.2, 5 % of 0.2–0.5, 2 % above 0.5 km². So the
  DEM cannot be the sole anchor; as a merge rule it only has to be right where two bodies share a depression.
- Known misses with this recipe on 19_39/2019: 26 of 217 Dunmire lakes untouched — 5 with no pixel ever > 0.5 and
  21 whose > 0.5 core is < 0.05 km² at 10 m (half are lakes < 0.1 km²). NDWI 0.3 recovers most but adds ~200
  slush components per tile; SR-based 0.5 recovers more but SR does not exist before 2019 over Greenland. Open:
  a 0.4 band, and site-level persistence across years (a small lake accumulates over ten seasons, a crevasse field
  moves).
- One polygon per component per year. Attributes: year, area, n_scenes_water, first/last water date.
- Domain: the whole GIMP ice mask, no elevation cap (lakes advance inland ~10 m a⁻¹, Fan 2025; the
  < 2000 m zone was only used to size the job). Minimum size is a definition of *lake*, not of water
  (Ryan 2026: features < 0.015 km² hold 38–67 % of surface water); lowering it later adds sites via the
  append rule without renumbering.

**Site.** The union of all outlines over the basis years, closed with a small kernel to bridge ice-lid gaps
(default closing radius 150 m; test 100–300 m on CW), then split into connected components.
A site is the "general region" a lake occupies. Sites do not change geometry within a registry version.
Precedent: Selmes et al. 2011 (GRL 38, L15501) summed five years of MODIS classifications into the maximum extent
of every lake and took one centre point per lake "from all five years of data, as some lakes did not form every
year", then tracked each lake inside that fixed site; they also tested (2000 vs 2010) that lakes do not advect.
We keep the union-as-site logic and the fixed centre as the label, and replace their seeded region-growing
(right at 250 m) with membership by overlap against the union polygon (right at 10 m).

**Lake ID.** Josh's convention (2026-09-05): `G00000_N691234_W0485678` = `G` for Greenland, a five-digit serial for
quick referral, then the registration centroid as hemisphere letter + fixed-width digits with four implied decimals
(69.1234 N, 048.5678 W). No signs or points. The serial is assigned once and never reused or reassigned: the first registry
numbers lakes north to south (within a drainage-basin sweep at the Greenland-wide build); each later batch of additions
(a new season, or a Landsat back-fill year) is appended after the last existing serial, ordered north to south within the
batch. The serial is for referral only — position is read from the coordinate part, never from the number. The coordinate part is a label frozen at registration, never recomputed
when a polygon is redrawn. (Tile-level test registries carry provisional serials.)

**Membership.** An outline belongs to every site it overlaps. Default: any overlap; report the fraction.
An outline overlapping ≥ 2 sites is an observation on each, flagged `shared`, area attributed by the part
inside each site with the remainder to the nearest site.

**Observation.** The site polygon is the basin delineation (largest extent water has ever revealed). Per site
and per scene date, from a sat-tile-stack cut around the site: (a) water pixels *inside* the site — threshold
0.18 as Dunmire 2025 used inside known outlines, since the outline already says "lake" and the question is only
how full; count, fraction of site, first/last date; (b) the usability flag of the scene window (cloudy-tile);
(c) *spill*: the total area of any water body overlapping the site that extends beyond it, so a frozen outline never
truncates a record (persistent spill → widen the site in the next versioned rebuild); (d) a body spanning ≥ 2 sites
gives each site its inside part and a `shared` flag. Per year, the seasonal-maximum outline row as in §1.

**Event.** `shared` observations are the coalescence record; drainage classification is a later layer with
its own method version.

## 2. Basis and versions

- v0.0 basis: Sentinel-2 seasons 2016–2025 (S2A from June 2015, S2B from March 2017; 2016 is thin).
  Dunmire's 2018 and 2019 outlines are the check for those two years.
- Registry v0.0 = all sites from the basis, frozen, with IDs. Released with a DOI and a crosswalk to
  Dunmire's `new_id` (per-year IDs → site ID) and to `CW_lookup.csv`.
- Append rule (new season, or Landsat back-fill of 1985–2015 applied backwards): build that year's outlines
  with the same recipe; each outline overlapping ≥ 1 site → observation(s); overlapping none → a new site
  registered with a new ID and `first_seen = year`. Sites are never merged or renumbered by the append rule.
- A rebuild with a longer basis is a new version (v1.0) with an old→new crosswalk; IDs carry over wherever
  the new site overlaps exactly one old site.

## 3. The DEM: an attribute, not a rule (Josh, 2026-09-05 evening)

The DEM is not a reliable anchor on its own: 9 % of 2018 and 12 % of 2019 Dunmire lakes in CW lie entirely outside any
32 m depression, 27 % of those under 0.1 km². So v0.0 does **not** merge or cut sites with the DEM. Sites come from the
ten-season union and the closing alone (§1), which already joins the crescent CW_0381 (`out/08_multiyear_sites_19_39.png`).
The DEM is computed *per site* and reported, never enforced: fraction of the site inside a depression, the sub-basin ID(s)
it overlaps, max fill depth, and whether the site's centroid falls in a depression. The paper then shows the correlation
(expected ~90 % of sites, ~98 % of sites > 0.5 km², in a depression) as evidence that sites are bed-controlled basins,
without forcing it. A DEM merge stays possible later as a versioned rebuild with a crosswalk.

Consequences: the one geometric knob is the site closing radius (default 150 m; sweep on the ten-season union before
freezing). The validation-by-labelling test (true-colour chips, ~100–200 cases) now targets the *closing* decisions —
pairs of per-year outlines that the closing joins, and near neighbours it leaves apart — rather than DEM merges.

## 3b. Deliberately out of v0.0 (additive later)

- Volume. Flagged by Josh 2026-09-05. Two hooks the basin-anchored design leaves ready: (i) an area→volume
  curve per site from the DEM hypsometry (lower bound from the mosaic, whose floor is a lid median; true floors
  from dated 2 m ArcticDEM strips after a drainage); (ii) radiative-transfer depths (Sneed & Hamilton 2007, used
  by everyone) need open water, so lids break them exactly where (i) helps; ICESat-2 crossings as validation
  (Datta & Wouters 2021; Melling 2025).
- Rivers, connectivity, drainage type, Sentinel-1 buried water.

## 4. Compute

- Finding sites is a per-pixel reduction over every scene of every season for every granule over the lake
  zone — the Earth Engine pattern (Dunmire 2021 and Fan 2025 ran exactly this there), not a sat-tile-stack
  job (that tool cuts stacks around *known* centroids). Order of magnitude: 100–150 granules × 40–100 scenes
  per season × 2 bands at 10 m ≈ 10–20 TB of band reads for ten seasons; server-side on GEE, nothing
  transferred. Exported annual 10 m masks (or polygons) are the reproducible artifact, kept with the code.
- Observing sites (per-date water fraction, later drainage type) is the sat-tile-stack job on the Planetary
  Computer, once sites exist.
- Sites, IDs, membership, crosswalks: seconds to minutes on the laptop once annual outlines exist.
- Landsat back-fill: Fan 2025's annual boundaries (TPDC, registration) if obtainable, else the same recipe
  on Landsat on GEE.

## 5. Tests before freezing (all light, CW, with what is on disk)

1. Closing radius sweep 0–300 m on the union of Dunmire 2018+2019: sites count, lid pairs joined, neighbours
   merged, versus the DEM sub-basin truth from `out/05_lake_to_basin_CW_p3.csv`.
2. Centroid-tolerance sweep 250 m–2 km on `out/pairs_2019_vs_2018.csv` (for the record: why overlap, not
   distance).
3. One granule, one season (2019, CW) on GEE with the recipe above, at k = 1, 2, 3 and with/without the
   opening, compared against Dunmire's 2019 outlines as a sanity check (not a target): pins k, the closing
   and the opening before the ten-season run.

## 6. Ten-season result on tile 19_39 (2026-09-05, `out/09_sites_19_39.*`)

Recipe as in §1 with the scene QA (§1, scene area > 3 × median of the top-10 scene areas dropped; it removed
the one striped 2024-08-31 acquisition that had painted 89 km² of diagonal bands into 2024 and nothing else).
Per-season outlines 167–207 (67–109 km²; 2019 the big year at 108.6 km²). Union + 150 m closing → **297 sites,
212 km²**. Versus Dunmire: 2018 139/139 touched, 138 covered ≥ 50 %; 2019 215/217 touched, 214 covered ≥ 50 %;
no Dunmire lake spans two sites; 7 (2018) and 12 (2019) sites hold two Dunmire lakes of the same year (the
crescent CW_0381 and dumbbell CW2018_1270 among them — one site each, from the union and closing alone).
Persistence is bimodal: 71 sites held water in all ten seasons, 53 in only one; small sites are the one-offs
(median 1 season below 0.1 km², 9–10 seasons above 0.5 km²). DEM correlation, reported not enforced: 66 % of
centroids fall in a 32 m depression overall, 36 % below 0.1 km², 89 % at 0.5–1 km², 85 % above 1 km², 100 %
above 5 km²; 24 % of sites lie entirely outside any depression. Registry IDs are serial north-to-south
(`G00000_…` at the top of the tile); polygons in `09_sites_19_39.geojson`, per-year presence in
`09_site_years_19_39.csv`, Dunmire crosswalk in `09_dunmire_crosswalk_19_39.csv`.
Open before freezing: the closing-radius labelling test (§3), the extras (sites touching no Dunmire lake in
2018/2019 — most are one-season small sites; some may be slush), and a second tile to check the recipe is not
tuned to 19_39.
