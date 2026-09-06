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
- Domain: the whole ice sheet, implemented as the BedMachine Greenland v6 150 m mask (grounded + floating ice; Morlighem et al. 2017, derived from the GIMP ice mask), nearest-upsampled to 10 m and applied to the water pixels before component labelling (`scripts/08g_ice_mask.py`, 2026-09-05; tile 19_39 is 98.8 % ice and the mask changes nothing there but one ice-marginal site's off-ice lobe); no elevation cap (lakes advance inland ~10 m a⁻¹, Fan 2025; the
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
208 km²** (with the ice-sheet domain mask). Versus Dunmire: 2018 139/139 touched, 138 covered ≥ 50 %; 2019 215/217 touched, 214 covered ≥ 50 %;
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

## 7. Second tile: 29_45 (NE Greenland, submitted 2026-09-05 20:30 PDT)

To check the recipe is not tuned to 19_39, the same ten seasons are being exported for tile 29_45 (78.0–79.2° N,
20.5–26.6° W: the NEGIS / Storstrømmen area, the ground of Hochreuther 2021, Lutz 2025 and Fan 2025, with a nunatak
zone and an ice-free coastal strip — 92.9 % ice by BedMachine, so the first tile where the domain mask matters).
Dunmire has 217 lakes there in 2019 (176 km²). Sentinel-2 orbit overlap at 79° N gives ~1 700–1 900 scenes per
season against ~640 at 70° N, so each export takes roughly three times longer. Driver: `scripts/08f_batch_tile.py`
(submits 08a for each season, polls, downloads with 08b); `scripts/chain_29_45.sh` runs 08e when the batch finishes;
log `out/08_batch_29_45.log`. Size threshold stays 0.05 km² (Dunmire 2021/2025, How 2025, Miles 2017's 0.0495 km²);
the NEGIS papers used 0.015 km², so a recall gap against them is expected and is a threshold choice, not a recipe fault.

### 1b. Two artefact rules added 2026-09-05 late (found by the closing test chips)

- **Swath-edge lines.** Three one-season "outlines" on tile 19_39 (2020 ×2, 2023 ×1, each ≤ 0.1 km²) were perfectly
  straight 1–2 px diagonal lines along Sentinel-2 swath edges (band misregistration at the footprint edge makes NDWI
  spike); a fourth was glued onto a real lake (G00277). Two rules: (i) at the source, 08a drops `EDGE_PX` = 3 pixels
  inside every scene's valid footprint (`focalMin`), for all future exports; (ii) in 08e/11a, a per-season outline must
  contain a 50 m wide core (survive erosion by a 2 px disk). On 19_39 the core rule removes exactly the three lines and
  touches no Dunmire lake; the attached line on G00277 needs the source rule (re-export pending).
- **Chip fetch honours scene QA.** 11b now excludes the scenes the tile's task JSONs dropped (the striped 2024-08-31
  scene had been picked as "wettest" for one chip).

## 8. Closing-accuracy test, tile 19_39 (2026-09-05 late; Josh and Claude labelled all 105 cases independently)

Cases (`out/11_cases_19_39.csv`; chips `out/11_chips/`; labels `out/11_labels_both_19_39.csv`; label page
https://claude.ai/code/artifact/b4b6d918-2fe0-4e47-a742-f5b1ad339d4c): 16 sites the 150 m closing glued from separate
union pieces ("joined"), 35 sites that held ≥ 2 outlines in one season ("lid", joined by the union, not the closing),
54 pairs of sites within 1 km left apart ("apart"). Answers 1 = one lake basin, 2 = two or more, 3 = can't tell.
Agreement 84/105 (80 %); disagreements are channel-linked ponds, slush fields, and the line artefacts (§1b).

| pool | n | Josh: wrong | Claude: wrong | reading |
|---|---|---|---|---|
| joined (closing acted) | 16 | 6 two-or-more | 10 two-or-more | the closing is wrong about half the time it acts |
| lid (union joined) | 35 | 4 | 3 | the union rule is right ~90 % |
| apart (left separate) | 54 | 6 "one lake" | 1 "one lake" | few false splits; Josh's six are at 155–951 m, no radius fixes them |

The bridged gap decides it (`out/11_joined_gaps_19_39.csv`): the five glued sites with gaps ≤ 80 m are lid remnants
and both labellers call them one lake (Josh 4/5, Claude 3/5 with one slush "can't tell"); of the eleven with gaps
114–201 m, Josh calls 6 and Claude 9 separate lakes. **Recommendation: site closing 50 m (5 px, the same radius as the
per-scene closing) instead of 150 m.** It keeps every lid join both labellers endorse, drops most of the wrong ones,
and creates no new false split at the gaps where the labellers disagree. Not applied yet (Josh's call; one 08e rerun
with `SITE_CLOSE_PX=5`). Registry after the core rule (§1b): 294 sites, Dunmire recall unchanged (139/139, 215/217).

### 8b. Decision, 2026-09-05 22:00 PDT (Josh): the appendage rule replaces the flat 150 m closing

Two bodies of water in the ten-season union are one site if (1) they are within 50 m of each other (the ice-lid case,
same radius as the per-scene closing), or (2) they are within 150 m and the smaller is less than one fifth of the
larger (an appendage: a tail, a spur, a pond hanging off a lake, which does not merit its own ID). Otherwise they are
separate sites with separate IDs. Inputs: gap and size ratio only. On tile 19_39 the three appendage joins are at 6–12 %
of their lake and the smallest genuine pair is at 30 %, so one fifth sorts every labelled case the way both labellers
did except the two Josh-one/Claude-two judgment calls. To implement in 08e (replace the single `SITE_CLOSE_PX`
closing): closing at 5 px → provisional sites; then for each pair of provisional sites within 15 px, merge if
min(area)/max(area) < 0.2; rebuild 19_39 and score against `out/11_labels_both_19_39.csv`; 29_45 inherits it.

### 8c. Appendage rule implemented and scored (2026-09-06, `scripts/08e` + `scripts/11e_score_rule.py`, `out/11_rule_scores_19_39.txt`)

08e no longer closes at the site level: union pieces are joined by distance (`JOIN_ALL_M`, default 50) or by distance and
size ratio (`JOIN_APP_M`, `APP_RATIO`); the joins are logged in `out/09_joins_{TILE}.csv`. Scored on the 105 labelled cases
(label 3 excluded per labeller; "agree" = the rule joins exactly the cases that labeller called one lake):

| setting | Josh | Claude | note |
|---|---|---|---|
| 50 m always, 150 m if < 1/5 | 88/103 (85 %) | 86/94 (91 %) | only 1 appendage join fires; two of Josh's three appendages sit at 187 and 190 m |
| 50 m, 200 m if < 1/5 | 90/103 (87 %) | 86/94 (91 %) | picks up both |
| 50 m, 250 m if < 1/5 | 91/103 (88 %) | 87/94 (93 %) | also joins the drained-lake moat pair G00258/G00256 that both labellers called one lake |
| old flat 150 m closing | 88/103 (85 %) | 83/94 (88 %) | for reference |

The ratio (1/5 vs 1/4) changes nothing on this tile. The remaining disagreements are pairs Josh joins at 300–950 m
(no distance rule reaches them) and three channel-linked lakes joined by the union itself (G00077, G00092, G00199),
which no distance rule can split. Registry on disk at the time of writing: 150 m setting, 306 sites, Dunmire recall
139/139 and 215/217 (one lake per year now spans two sites). Awaiting Josh's pick of the appendage distance.

## 9. Tile 29_45 result (NE Greenland, 2026-09-06, built with the old 150 m closing + core rule; `out/09_sites_29_45.*`)

Ten seasons, 1 643–2 409 scenes each; per-season area swings from 38 km² (2018) to 185 km² (2023). 309 sites, 359 km².
Dunmire recall: 2018 108/109, 2019 217/217, none split. Persistence bimodal again (52 sites all ten seasons, 40 one).
DEM correlation weaker than CW: 55 % of centroids in a 32 m depression (66 % on 19_39), 34 % of sites entirely outside
one (24 %); by size 20 % below 0.1 km² rising to 76 % above 1 km². Elongated NE–SW sites along the NEGIS shear margin
are real flow-stripe lakes. The recipe is not tuned to 19_39: recall is as good on a tile with a different climate,
orbit geometry and terrain. Rebuild with the appendage rule once its distance is chosen.

## 10. Visual check and two review flags (2026-09-06 afternoon; `scripts/12a_truecolor_fetch.py`, `12b_truecolor_figs.py`)

Josh's decision: appendage rule at 50 m / 250 m / one fifth (08e defaults). Both registries rebuilt: 19_39 303 sites
(5 lid + 4 appendage joins), 29_45 320 sites (8 + 11); Dunmire recall unchanged. Visual check: six 12 km windows per
tile where sites are densest, each on the wettest clear scene of its two wettest seasons (cloud probability inside the
window < 20 %, July–mid-September, water counted only inside the window's own outlines), site outlines in red with the
serial on every polygon part (`out/12_truecolor/{TILE}_w{k}.png`). Reading: on clear windows every outline is one lake
and the outline holds while the water inside changes between years; ten-season outlines include slush aprons around
lakes in slush zones (jagged perimeters, e.g. 19_39 serials 238, 247, 283) — the water-seen footprint, by design.
Two per-site review flags, attributes not rules: `ice_marginal` (site within 300 m of non-ice on the BedMachine mask;
2 on 19_39, 19 on 29_45 — fjord and ice-dammed water at nunataks and the coast, whether they belong in a supraglacial
registry is a scope call), and `thin` (widest point ≤ 60 m; catches line-like sites such as 19_39 serial 214, a
100 m wide, 1.5 km straight strip that passed the 50 m core rule). `elongation` (P²/4πA) is reported too but measures
jaggedness, not linearity.

## 11. Two more rules decided 2026-09-06 (Josh): minimum width 40 m; sites touching rock or ocean are out

- **Minimum width, per-season outline** (`MINW_PX` = 2 → 40 m in 08e): open the outline with a 2 px disk, then grow the
  survivors back by 2 px inside the original outline. Lake shapes are unchanged; channels narrower than 40 m and pond
  chains necked below 40 m are cut; pieces under 0.05 km² drop. Tested (`scripts/13_min_width_test.py`,
  `out/13_min_width_*`): at 40 m no Dunmire lake is lost on either tile and ~0.5 % of water goes; at 60 m the NE tile
  loses four Dunmire lakes that are 50–60 m wide flow-stripe troughs, so 40 m was chosen. Aspect ratio was rejected:
  a pond chain is not elongated as a whole and a winding channel can have any ratio; width is what both lack.
- **Ice-marginal exclusion** (`EXCLUDE_TOUCHING=1`): a site whose water touches the non-ice classes of the BedMachine
  mask (ice-free land or ocean) is an ice-marginal lake (How et al. 2025's domain) and is excluded from the registry;
  sites within 300 m of non-ice that do not touch it stay, flagged `ice_marginal`. Before the rule: 1 touching site on
  19_39, 11 on 29_45 (20.5 km², the 9.7 km² fjord body among them).
Result on 19_39 with both rules: 301 sites (5 lid + 3 appendage joins), 194 km²; Dunmire 139/139 and 215/217 unchanged;
the 100 m × 1.5 km strip (old serial 214) is gone and no site is flagged thin; the one touching site (0.1 km²) excluded.
**Reversed the same evening (Josh):** the 40 m rule trims thin parts *attached* to lakes (outflow tails, narrow moat
segments), which are part of the lake object. What Josh wants is "drop an outline only if it is thin everywhere; keep
whole anything that has a wide body" — and that is the existing 50 m core rule (§1b): a component survives, complete,
if any part of it is ≥ 50 m wide, so a lone channel goes and a channel attached to a lake stays with the lake. `MINW_PX`
now defaults to 0 (knob kept); the core rule is the width rule. Pond chains with wide ponds remain one object.

## 12. Compute route for Greenland-wide (decided 2026-09-06 evening, Josh): Level-1C from Copernicus Data Space, on Sherlock

Earth Engine's noncommercial tier (150 EECU-hours/month; two tiles × ten seasons spent it) cannot carry the ~67 lake-zone
tiles. Decision: the counts (08a's reduction) run on Sherlock as one array task per tile-season, reading Sentinel-2
**Level-1C** (Dunmire's product, and the one our Earth Engine exports use) from the Copernicus Data Space `eodata` S3
store (free account; JPEG2000 in SAFE, so whole-band reads). Level-2A was rejected for the counts: (i) atmospheric
correction changes the blue band most, so the 0.5 NDWI threshold tuned on top-of-atmosphere does not carry over;
(ii) it is incomplete before 2019 (July over CW: Planetary Computer 13/24/25 scenes in 2016/17/18, AWS 0/20/34, against
651/572/627 Level-1C scenes per whole season on the tile), and the basis needs every season. Level-2A stays fine for
2019+ products (observation layer, chips). Planetary Computer serves Level-2A only; AWS Level-1C is requester-pays.
Test plan before scaling (none started; needs Josh's Copernicus account keys on Sherlock, never in the repo):
1. one scene on one node (auth, decode, throughput); 2. tile 19_39, season 2019, one node, output in 08b's format,
site builder unchanged, compare outlines and Dunmire recall with the Earth Engine result, record wall-clock and
bytes; 3. ten tile-seasons as a ten-node array job (throttling, Greenland-scale estimate). No 20 m or per-day cuts:
parallelism covers the volume.
