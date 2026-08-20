/**** =====================================================================
 *  SOURCE: https://code.earthengine.google.com/e32121007cd578cd8c741d9a7815c02f?noload=1
 *  Referenced in: 4_Typology_OLDARCHITECTURE/1_makeDyad.R:22
 *  Retrieved:     2026-08-12
 *  Produces:      outCntr_Typol.csv, outTBA_Typol.csv,
 *                 outB{5,10,100,200}k_Typol.csv  -> Drive folder TBA2
 *                 consumed by 4_typology/1_data/geeOut/
 *
 *  IMPORTS (9). The Code Editor stores imports with the link, not in the
 *  source text. Recorded here so this file is self-describing. Do NOT paste
 *  these into a script that already has the imports panel populated.
 *
 *    var Irr1k   = ee.Image('projects/globaltbav2/assets/IrrigationProjections/MaierIrrigatedArea');
 *    var GWf     = ee.Image('projects/globaltbav2/assets/IrrigationProjections/GWf');
 *    var GWS     = ee.Image('projects/globaltbav2/assets/IrrigationProjections/GWS_Rosa2020');
 *    var BWS     = ee.Image('projects/globaltbav2/assets/IrrigationProjections/BWS_Rosa2020');
 *    var urb     = ee.ImageCollection('JRC/GHSL/P2023A/GHS_SMOD');
 *    var Crop    = ee.Image('projects/globaltbav2/assets/IrrigationProjections/Crop1k_0307');
 *    var TBA     = ee.FeatureCollection('projects/globaltbav2/assets/VectorBuffers/_TBA_Split_2025');
 *    var lb      = ee.FeatureCollection('projects/globaltbav2/assets/VectorBuffers/_LandBorders_2025');
 *    var Country = ee.FeatureCollection('projects/globaltbav2/assets/VectorBuffers/GADM_Simplified');
 *
 *  KNOWN ISSUE (see Notes/typology_audit.md and gee/dyad_border_corrections.js):
 *  the four distance-buffer exports measure distance to `lb` — ANY land
 *  border — and reduce over TBA segments. A segment therefore carries one
 *  near-zone area reused across every dyad it enters. Fixed downstream by
 *  dyad_border_corrections.js, not here.
 *
 *  DOC/CODE MISMATCH: the header block below says CR3/IR3/GW3 are all
 *  "groundwater stress", and the variable names all read *_GWS3, but the
 *  code masks CR3 on GWS and IR3/GW3 on BWS. The inline comments are the
 *  accurate ones. Names and header are misleading; the SI's Overdraft axis
 *  is defined on blue-water scarcity, so the BWS masking looks deliberate.
 * ==================================================================== */

/**** ---------------------------------------------------------------------
Groundwater-dependent irrigation and cropland statistics
--------------------------------------------------------------------------
This script aggregates several spatial indicators at administrative or
basin levels using Google Earth Engine.
Pixel-level variables are converted to area (hectares × 10³ = kHa) using
pixel area and then summed over polygons using reduceRegions.
Main indicators produced:
- CR   : Cropland area
- IR   : Irrigated cropland area
- GW   : Groundwater-irrigated cropland
- UR   : Urban area
- CR3  : Cropland where groundwater stress > threshold
- IR3  : Irrigated cropland where groundwater stress > threshold
- GW3  : Groundwater irrigation where groundwater stress > threshold
Outputs are summarized for:
1. Countries
2. TBA regions
3. TBA regions within distance buffers of land boundaries
   (5 km, 10 km, 100 km, 200 km)
All outputs are exported as tables to Google Drive.
---------------------------------------------------------------------- ****/
/**** ---------------------------------------------------------------------
Prepare base rasters
---------------------------------------------------------------------- ****/
// Ensure irrigation mask has no missing values
Irr1k = Irr1k.gt(0).unmask()
// Cropland area in kHa per pixel
var CR = Crop.unmask()
  .multiply(ee.Image.pixelArea())
  .divide(1e7)
// Tile scale used in reducers (controls memory usage)
var tiles = 16
/**** ---------------------------------------------------------------------
Urban area
---------------------------------------------------------------------- ****/
// Urban area derived from urban dataset (year 2005)
// Threshold >20 represents urban class
var Urb = ee.Image(
    urb.filter(ee.Filter.calendarRange(2005,2005,'year')).first()
  )
  .gt(20)
  .multiply(ee.Image.pixelArea())
  .divide(1e7)
/**** ---------------------------------------------------------------------
Derived irrigation and groundwater indicators
---------------------------------------------------------------------- ****/
// Cropland located in groundwater stress areas
var CR_GWS3 = CR.multiply(GWS.gt(2))
// Irrigated cropland
var IR = Irr1k.unmask().multiply(CR)
// Groundwater-irrigated cropland
var GW = Irr1k.multiply(GWf).multiply(CR)
// Groundwater irrigation where basin water stress > threshold
var GW_GWS3 = GW.multiply(BWS.gt(2))
// Irrigated cropland where basin water stress > threshold
var IR_GWS3 = IR.multiply(BWS.gt(2))
/**** ---------------------------------------------------------------------
Create image stack used for reductions
---------------------------------------------------------------------- ****/
var IM = ee.Image.pixelArea().rename('Area')
  .addBands(GW.rename('GW'))
  .addBands(CR_GWS3.rename('CR3'))
  .addBands(Urb.rename('UR'))
  .addBands(GW_GWS3.rename('GW3'))
  .addBands(IR_GWS3.rename('IR3'))
  .addBands(CR.rename('CR'))
  .addBands(IR.rename('IR'))
/**** ---------------------------------------------------------------------
Country-level aggregation
---------------------------------------------------------------------- ****/
var outCntr = IM.reduceRegions({
  collection: Country,
  reducer: ee.Reducer.sum(),
  scale: 1000,
  tileScale: tiles
})
.map(function(ft){ return ft.setGeometry(null) })
Export.table.toDrive(
  outCntr,
  'outCntr_Typol',
  'TBA2',
  'outCntr_Typol'
)
/**** ---------------------------------------------------------------------
TBA regional aggregation
---------------------------------------------------------------------- ****/
Map.addLayer(TBA)
// Sum statistics over TBA polygons
var outTBA = IM.reduceRegions({
  collection: TBA,
  reducer: ee.Reducer.sum(),
  scale: 1000,
  tileScale: tiles
})
.map(function(ft){ return ft.setGeometry(null) })
Export.table.toDrive(
  outTBA,
  'outTBA_Typol',
  'TBA2',
  'outTBA_Typol'
)
/**** ---------------------------------------------------------------------
Distance-to-boundary analysis
Pixels are limited to areas within specified distances from lb
---------------------------------------------------------------------- ****/
// 10 km buffer
var outB10k = IM
  .multiply(lb.distance(10000).lt(10000))
  .reduceRegions({
    collection: TBA,
    reducer: ee.Reducer.sum(),
    scale: 1000,
    tileScale: tiles
  })
  .map(function(ft){ return ft.setGeometry(null) })
Export.table.toDrive(outB10k,'outB10k_Typol','TBA2','outB10k_Typol')
// 5 km buffer
var outB5k = IM
  .multiply(lb.distance(5000).lt(5000))
  .reduceRegions({
    collection: TBA,
    reducer: ee.Reducer.sum(),
    scale: 1000,
    tileScale: tiles
  })
  .map(function(ft){ return ft.setGeometry(null) })
Export.table.toDrive(outB5k,'outB5k_Typol','TBA2','outB5k_Typol')
// 100 km buffer
var outB100k = IM
  .multiply(lb.distance(100000).lt(100000))
  .reduceRegions({
    collection: TBA,
    reducer: ee.Reducer.sum(),
    scale: 1000,
    tileScale: tiles
  })
  .map(function(ft){ return ft.setGeometry(null) })
Export.table.toDrive(outB100k,'outB100k_Typol','TBA2','outB100k_Typol')
// 200 km buffer
var outB200k = IM
  .multiply(lb.distance(200000).lt(200000))
  .reduceRegions({
    collection: TBA,
    reducer: ee.Reducer.sum(),
    scale: 1000,
    tileScale: tiles
  })
  .map(function(ft){ return ft.setGeometry(null) })
Export.table.toDrive(outB200k,'outB200k_Typol','TBA2','outB200k_Typol')
