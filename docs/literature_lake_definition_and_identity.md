# How the literature defines "a lake" and "the same lake in another year"

*Compiled 2026-09-05 (Pacific) from the PDFs in Josh's Zotero library plus open-access fetches; every
claim below was read in the paper's own methods text, not taken from memory. Sources not fetched in
full are marked (abstract only). Complements README §1 (Dunmire's recipe) and §3 (prior art).*

## A. What a lake is (detection layer)

| Paper | Sensor / grid | Water rule | Minimum object | Season aggregation |
|---|---|---|---|---|
| Sundal 2009 RSE | MODIS 250 m | band thresholds (Box & Ski style) | ~0.1 km² | per image |
| Selmes 2011 GRL, 2013 TC | MODIS 250 m | red < 65 % of surrounding ice (dynamic) | 0.125 km² = 2 px, reached at least once in 5 yr | daily areas per lake; annual max per lake |
| Liang 2012 RSE | MODIS 250 m | histogram-slope dynamic threshold on reflectance difference in a 25×25 window | 2 px, seen ≥ 3 times | daily; within-season matching |
| Morriss 2013 TC | ETM+ 30 m + MODIS 250 m | NDLI (red−NIR)/(red+NIR) > 0.21 | 0.125 km² and 10-yr max depth > 2 m | 78 fixed basins |
| Fitzpatrick 2014 TC | MODIS 250 m | modified NDWI + red/blue band limits, > 50 % water pixels | 0.0625 km² = 1 px | daily; per-year max |
| Arnold 2014 TC / Banwell 2014 | Landsat | blue/red ratio > 3 (Paakitsoq: 1.4 in Banwell) | – | per image |
| Williamson 2017 RSE (FAST) | MODIS 250 m | red band MOD09 > 0.640 dynamic; water on ≥ 3 occasions | 0.125 km² at least once | seasonal maximum-extent polygons; areas tracked inside them |
| Miles 2017 Front. | Landsat 8 30 m | NDWI (blue−red)/(blue+red) > 0.5 (empirical) | 55 px = 0.0495 km² at least once | composite max mask dilated 2 px; FAST inside |
| Cooley & Christoffersen 2017 JGR | MODIS 250 m | spatially averaged normalized reflectance difference | > 0.25 km²; ≥ 3 images, twice in 6 d | yearly water mask → potential lakes; present if > 20 % or > 2 px wet |
| Williamson 2018 TC | S2 10 m + L8 30 m | NDWI + band thresholds (Moussavi/Pope style) | < 0.125 vs ≥ 0.125 km² classes | seasonal max extents; FAST |
| Gledhill & Williamson 2018 Ann. Glac. | Landsat 30 m, 1985–2016 | manual delineation | 0.0009 km² (1 px) | July–Aug max extent per year |
| Moussavi 2020 RS (Antarctic) | L8 / S2 | NDWI(blue,red) > 0.19 (L8) / 0.18 (S2) + green−red, blue−green limits | < 5 px and < 2 px wide removed | per image |
| Yang 2021 J. Glac. | L8 30 m | NDWI_ice > 0.25 global; width > 5 px = lake, else river | – | per image |
| Dunmire 2021 TC / 2025 ESS | S2 at 30 m | NDWI(blue,red) > 0.5; SWIR > 0.140 cloud; NDSI < 0.85 & blue < 0.4 margin; two 5×5 closings | > 0.05 km² | one seasonal-max outline per lake-year; 2025 adds p_water inside it (NDWI > 0.18) |
| Turton 2021 TC (79N) | S2 10 m | blue/red ratio > 1.6 | 150 px = 0.015 km² | per date |
| Hochreuther 2021 RS (NEGIS) | S2 10 m | blue/red ratio > 1.6, shadow model, sieve | 150 px = 0.015 km² | per date, 1.5 d spacing 2016–2019 |
| Hu 2022 RS (abstract only) | S2 | random forest with texture/morphology | – | per melt season 2016–2018 |
| Otto 2022 Front. (Ryder, 1985–2020) | Landsat 15–30 m | NDWI, image-specific threshold 0.15–0.5 | – | July–Aug masks merged per year |
| Lutz 2023 RS (abstract only) / Lutz 2025 TC | S2 10 m | U-Net segmentation | – | near-daily 2016–2022 |
| Wang & Sugiyama 2024 RSE | S2 + L8 at 10 m | random forest on bands + indices | ≤ 18 px clusters and ≤ 3 px wide removed | per image; 8-yr cumulative max |
| Fan 2025 J. Glac. (1985–2023) | Landsat 30 m | per-year max-NDWI_ice composite (qualityMosaic), NDWI_ice > 0.2 | < 5 px (4500 m²) and < 2 px wide removed; manual correction | one July–Aug max extent per year |
| How 2025 ESSD (ice-marginal) | S1 + S2 + ArcticDEM sinks | three methods, dissolved | ≥ 0.05 km² | per inventory year |
| Dean 2026 TC (NE Greenland winter drainage) | L8/9 30 m + S1 | NDWI_ice > 0.4 | ≥ 0.1 km² | 10-yr composite of summer outlines |
| Ryan 2026 AGU Adv. (SkySat) | 0.5 m + S2 | deep learning; NDWI > 0.10 / 0.20 checks | none; shows < 0.015 km² features hold 38–67 % of water | per image |
| Stokes 2019 Sci. Rep. (EAIS) | L8 / S2 | NDWI(green,NIR) > 0.3 | 2 px | January 2017 only |
| Arthur 2022 Nat. Comm. (EAIS) | L8 30 m | NDWI, < 5 px and < 2 px wide removed | – | January max mask per year |
| GLAKES, Pi 2022 Nat. Comm. (global, non-glacial) | Landsat | water-occurrence probability per period | > 0.03 km² max extent | three periods; 1° cells |

Ignéczi 2016 GRL (DEM only): GIMP 30 m sink fill; drop ≤ 0.125 km², mean depth ≤ 1.5 m or ≥ 50 m, thin/floating ice;
recall 75 % of observed lakes below the ELA sit in a depression, precision 19 % of depressions host a lake.
Arnold 2014: 78 % of lake centroids fall in DEM depressions. Otto 2022: GIMP Cut-Fill sinks, "nearly all" 2020 lakes in one.

## B. What "the same lake" is across time (identity layer)

Four families, from weakest to strongest identity:

1. **None — per-year objects, per-year IDs.** Dunmire 2021/2025 (`CW2019_1524`; the 2025 paper compares
   "lakes that appear in the same topographical depression in both years" for case studies only), Cooley &
   Christoffersen 2017 ("each lake in each year is considered a unique observation"), Gledhill & Williamson 2018,
   Otto 2022 ("all individual lakes in each annual lake mask were delineated"), Hu 2022, Sundal 2009, Stokes 2019.
2. **Per-pixel recurrence — no lake objects across years.** Fan 2025 (reoccurrence = number of years a pixel is
   in the max composite, 1985–2023), Wang & Sugiyama 2024 (recurrence frequency of water pixels 2014–2021),
   Dirscherl 2021 (January recurrence 2016–2021, notes ice-flow advection bias), Arthur 2022 ("Count Overlapping
   Features" on max-extent polygons), Qiu & Ran 2023/2025 (fixed occurrence grids), GLAKES (max-extent boundary +
   1° cells). Liang 2012 sits between 1 and 2: within a season lakes are matched by overlap and similar area
   (0.5·A1 < A2 ≤ 2·A1 + 8 px); across years a "recurrence frequency" (years active / 10) is reported per lake
   (13.2 % reappear every year) without stating the cross-year matching rule.
3. **Fixed sites from the multi-year water record.** Selmes 2011: sum 5 years of classifications → maximum extent
   of every lake → one centre point per lake ("centres were taken from all five years of data, as some lakes did not
   form every year") → each MODIS image classified by region-growing from those seeds; 2 600 lakes; tested no
   advection between 2000 and 2010. Morriss 2013: 10-yr max depth > 2 m and > 0.125 km² → 78 fixed lake basins
   with a 500 m buffer. Williamson 2017 FAST / Miles 2017 / Stevens 2026 FASTER: seasonal maximum-extent
   polygons are the units inside which areas are tracked (within one season). Dean 2026: all summer outlines of
   ten years merged into one composite of maximum extents → 404 fixed lake polygons, "so lakes that remain
   partially or fully buried over multiple years" are still counted. How 2025 (ice-marginal): IDs "initially
   defined as overlapping water bodies across all inventory years", then manual common IDs for lakes made of
   several polygons, dissolved by ID into a maximum extent. Fitzpatrick 2014: "21 % reoccurred annually",
   "most lakes form in the same locations each year" (rule not stated).
4. **DEM sinks as the lake unit, with IDs.** Hochreuther 2021: ArcticDEM gridded to 100 m, WhiteboxTools
   fill_depressions, filled − original < 0 = sink, 3×3 majority filter, "assigned an explicit ID to each sink";
   "to be able to track individual lakes, a spatial join was conducted to attribute the lakes to their respective
   sinks" → 1 035 sinks, 880 water-filled at least once in 2016–2019. A water polygon overlapping more than one
   sink is treated as *cloud* (false-positive filter), i.e. coalescence is explicitly not modelled. Lutz 2023
   (deep-learning time series over the same area) and Lutz 2025 inherit those IDs: "each lake in our area of
   interest was assigned an identification number. Tracking specific lakes is possible since the surface
   depressions where lakes form are directly influenced by bedrock topography and thus stay in the same location";
   152 lakes drained at least once in 2016–2022, drainage wheels per lake ID. Ignéczi 2016 and Leeson 2015 use DEM
   depressions as *potential* lakes (site inventory, no water assignment).

## C. Reading

- The physical premise everyone cites is Echelmeyer 1991 / Lampkin & Vanderberg 2011 / Sergienko 2013: lakes sit in
  bed-controlled surface depressions and do not advect. Selmes 2011 tested it (2000 vs 2010: no advection).
- Nobody publishes a Greenland-wide, multi-year lake ID set. The two regional precedents that do carry IDs across
  years are Hochreuther 2021 → Lutz 2025 (DEM sinks, NEGIS, 2016–2022) and Dean 2026 (10-yr composite polygons,
  NE Greenland); the ice-marginal inventory (How 2025) is the only one with a stated overlap-plus-curation ID rule.
- Josh's hybrid is family 3 + family 4 combined; both halves have precedent, the combination and the coalescence
  event do not. Hochreuther's "overlaps more than one sink = cloud" is the exact case the hybrid handles
  differently.
- Minimum-object thresholds in use: 0.0625 / 0.125 / 0.25 km² (MODIS pixel counts), 0.015 km² (S2, 150 px),
  0.05 km² (Dunmire, How), 0.1 km² (Dean). Ryan 2026 shows features < 0.015 km² hold 38–67 % of surface water
  in May/August, so any threshold is a choice about *lakes*, not about water.
- Season aggregation is nearly universal: one July–August (Greenland) or January (Antarctica) maximum extent per
  year for interannual work; daily/near-daily areas inside fixed units for drainage work.

## D. Files
Extracted texts and grep hits are in the session scratchpad only (not kept). PDFs: Zotero storage, keys in
`claudiary/20260905S.md`. Open-access fetches: Turton 2021, Dean 2026, How 2025, Otto 2022, Miles 2017,
Stokes 2019, Arthur 2022, Pi 2022 (web pages); Hochreuther 2021 (AWI EPIC PDF).
