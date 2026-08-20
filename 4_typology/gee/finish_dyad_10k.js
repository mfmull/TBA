/**** ---------------------------------------------------------------------
FINISH THE 10 km DYAD ZONES — the 6 aquifers that never landed
---------------------------------------------------------------------------
100 km is complete (594/594). 10 km has 500 of 594 rows; 94 are missing
across six aquifer codes. This script exports only those, only at 10 km.

WHY THE ORIGINAL TASK HUNG, PROPERLY THIS TIME
----------------------------------------------
Two causes stacked, and the second is the one that matters:

 1. toList(CHUNK, k*CHUNK) on a computed collection is O(offset), so the last
    chunk was always the slowest. Fixed by partitioning on aquifer code.

 2. The last chunk happened to hold S013 — the Amazon system, 3.7 million km²
    over six riparians, with the Brazilian segment alone at 1.8 million km².
    Its zone geometry is
        segment_BRA ∩ segment_COL.buffer(10 km)
    and because that geometry is COMPUTED, reduceRegions re-derives it for
    every raster tile it touches. An intersection between two continental
    multipolygons, recomputed per tile, does not finish.

So the fix is not smaller chunks. It is: do the vector work ONCE, write it to
an asset, and reduce over the materialised result. The zone itself is a thin
10 km ribbon — a few thousand km² — so the reduction is trivial once the
geometry is no longer being recomputed.

HOW TO RUN
----------
Just run it. USE_ASSETS = false exports straight to Drive, 11 tasks, one
pass. Five of these six aquifers are 16,000-41,000 km² and were never the
problem; the two-stage detour is not worth imposing on them.

The escape hatch is there because an EE asset is the only thing GEE can read
back as a static, indexed FeatureCollection — a Drive CSV returns geometry as
WKT in a .geo column that ee.FeatureCollection() cannot consume, so "write
geometry to Drive and reduce it" would need a manual re-upload and produce
the same asset anyway. If and only if the S013 tasks fail, set
USE_ASSETS = true and run STAGE 1 then STAGE 2 for those.

Files land as dyadB10k_r1_<CODE>.csv and go straight into
4_typology/1_data/geeOut/ alongside the parts already there.

NOTE ON THE 4 OVERLAPPING ROWS
------------------------------
EU095 - EU098 already contributed 4 rows (BLR-LTU, BLR-POL, BLR-RUS,
LTU-BLR) to the parts that did finish. Re-running the code re-exports them,
so those keys will appear twice. apply_dyad_zones() now de-duplicates keys
whose values agree and stops only where they disagree, so this is harmless —
do not hand-edit the CSVs.
---------------------------------------------------------------------- ****/

/**** ---------------------------------------------------------------------
Configuration
---------------------------------------------------------------------- ****/

// false: one pass, straight to Drive. This is the default and should work.
// true : two passes for the codes in ASSET_CODES only — run STAGE 1, wait,
//        then STAGE 2. Use only if a Drive task actually fails.
var USE_ASSETS = false;
var STAGE      = 1;
var ASSET_CODES = ['S013'];   // which codes take the two-stage route

var D   = 10000;            // 10 km only. 100 km is already complete.
var TAG = '10k';
var RUN_TAG = 'r1';         // keeps these files distinct from the earlier ones

// The six aquifers still missing rows at 10 km.
var CODES = ['EU095 - EU098','EU099','EU282','S013','S021','S022'];

// Codes exported one task PER COUNTRY rather than one per aquifer, so a
// single continental segment cannot take the whole aquifer down with it.
// S013 is the Amazon; its six segments total 3.7 million km².
var SIDE_SPLIT_CODES = ['S013'];

var ASSET_DIR    = 'projects/globaltbav2/assets/VectorBuffers';
var ASSET_PREFIX = 'dyadZones10k_';

var P_CODE = 'code';
var P_AQ   = 'aq_id';
var P_ISO  = 'CC';

var VARS  = ['Area','CR3','GW','GW3','IR','IR3','UR'];
var tiles = 8;              // more headroom than the 4 used elsewhere

var SIMPLIFY_M = 500;       // counterpart only; never the reduction region
var MAX_ERROR  = 100;

/**** ---------------------------------------------------------------------
Base rasters — identical construction to the production script
---------------------------------------------------------------------- ****/
Irr1k = Irr1k.gt(0).unmask()

var CR = Crop.unmask()
  .multiply(ee.Image.pixelArea())
  .divide(1e7)

var Urb = ee.Image(
    urb.filter(ee.Filter.calendarRange(2005,2005,'year')).first()
  )
  .gt(20)
  .multiply(ee.Image.pixelArea())
  .divide(1e7)

var CR_GWS3 = CR.multiply(GWS.gt(2))
var IR      = Irr1k.unmask().multiply(CR)
var GW      = Irr1k.multiply(GWf).multiply(CR)
var GW_GWS3 = GW.multiply(BWS.gt(2))
var IR_GWS3 = IR.multiply(BWS.gt(2))

var IM = ee.Image.pixelArea().rename('Area')
  .addBands(CR_GWS3.rename('CR3'))
  .addBands(GW.rename('GW'))
  .addBands(GW_GWS3.rename('GW3'))
  .addBands(IR.rename('IR'))
  .addBands(IR_GWS3.rename('IR3'))
  .addBands(Urb.rename('UR'))

/**** ---------------------------------------------------------------------
Helpers
---------------------------------------------------------------------- ****/
var nRipDict = ee.Dictionary(TBA.aggregate_histogram(P_CODE)).getInfo();

var withN = TBA.map(function(f){
  return f.set('__n', ee.Dictionary(nRipDict).getNumber(f.get(P_CODE)))
})

// Export descriptions and asset ids allow only [A-Za-z0-9_-], and the R-side
// pattern ^dyadB10k(_[A-Za-z0-9]+)*[.]csv$ wants single-underscore groups.
// 'EU095 - EU098' -> 'EU095EU098'.
function slug(s){ return s.replace(/[^A-Za-z0-9]/g, '') }

function unitsFor(code) {
  if (SIDE_SPLIT_CODES.indexOf(code) < 0) return [{code: code, cc: null}];
  var ccs = withN.filter(ee.Filter.eq(P_CODE, code))
                 .aggregate_array(P_ISO).distinct().sort().getInfo();
  return ccs.map(function(cc){ return {code: code, cc: cc} });
}

// 'a' is the segment being measured, 'b' the counterpart. Restricting the
// LEFT side before the join is what makes a per-country task cheap: the join
// never builds the rows this task is not exporting.
function sidesOf(u) {
  var all = withN.filter(ee.Filter.eq(P_CODE, u.code));
  var left = u.cc ? all.filter(ee.Filter.eq(P_ISO, u.cc)) : all;
  return ee.Join.inner('a','b').apply(left, all, ee.Filter.and(
    ee.Filter.equals({leftField: P_CODE,  rightField: P_CODE}),
    ee.Filter.notEquals({leftField: P_ISO, rightField: P_ISO})));
}

function nameOf(u) {
  return 'dyadB' + TAG + '_' + RUN_TAG + '_' + slug(u.code) +
         (u.cc ? '_' + u.cc : '');
}
function assetOf(u) {
  return ASSET_DIR + '/' + ASSET_PREFIX + RUN_TAG + '_' + slug(u.code) +
         (u.cc ? '_' + u.cc : '');
}

/**** ---------------------------------------------------------------------
Zone geometry

zone = segment_A ∩ segment_B.buffer(d)

A pixel in A within d of B is within d of the A-B border: reaching B from A
means crossing it, and the nearest point of B from inside A lies on the
shared boundary. Clipping to the aquifer is automatic — both segments already
are the aquifer.

The counterpart is first clipped to the bounding box of A grown by d. Any
point of B outside that box is further than d from every point of A, so the
clip is EXACT — it cannot change the answer — and it stops the buffer from
being taken over a continental polygon when only a border strip can matter.
---------------------------------------------------------------------- ****/
function zonesOf(u) {
  return sidesOf(u).map(function(row){
    var A  = ee.Feature(row.get('a'));
    var B  = ee.Feature(row.get('b'));
    var gA = A.geometry();
    var reach = gA.bounds(MAX_ERROR).buffer(D, MAX_ERROR);
    var gB = B.geometry().intersection(reach, MAX_ERROR)
                         .simplify(SIMPLIFY_M);
    var zone = gA.intersection(gB.buffer(D, MAX_ERROR), MAX_ERROR);
    return ee.Feature(zone, {
      code:     A.get(P_CODE),
      aq_id:    A.get(P_AQ),
      CC:       A.get(P_ISO),
      CC_other: B.get(P_ISO),
      n_rip:    A.get('__n'),
      zone_m2:  zone.area(MAX_ERROR)
    });
  });
}

// A zone is EMPTY whenever the two segments are more than D apart. That is
// common and it is not an error: 142 of the 594 pairs are more than 10 km
// apart, BLR-RUS on EU095 - EU098 by 65 km. Their correct 10 km value is
// ZERO, and the row must still be exported -- a dyad missing from the CSV
// keeps its lb-based segment value in apply_dyad_zones(), which is exactly
// the number being corrected. Dropping empties would silently un-correct
// nearly a quarter of the set.
//
// reduceRegions tolerates empty geometry (the finished parts carry 115
// zero-Area rows). Export.table.toAsset does not -- "Unable to export
// features with empty geometry" -- so the asset route has to split them out
// and write the zeros directly.
function splitZones(z) {
  var zeros = ee.Dictionary.fromLists(VARS,
    ee.List.repeat(0, VARS.length));
  return {
    full:  z.filter(ee.Filter.gt('zone_m2', 0)),
    empty: z.filter(ee.Filter.lte('zone_m2', 0))
            .map(function(f){ return f.setGeometry(null).set(zeros) })
  };
}

/**** ---------------------------------------------------------------------
Export
---------------------------------------------------------------------- ****/
var units = [];
CODES.forEach(function(c){ units = units.concat(unitsFor(c)) });

print('stage', STAGE, '  tasks', units.length,
      units.map(function(u){ return u.cc ? u.code + '/' + u.cc : u.code }));

function reduceToDrive(fc, u) {
  var out = IM.reduceRegions({
    collection: fc,
    reducer:    ee.Reducer.sum(),
    scale:      1000,
    tileScale:  tiles
  }).map(function(f){ return f.setGeometry(null) });
  Export.table.toDrive(out, nameOf(u), 'TBA2', nameOf(u));
}

units.forEach(function(u){
  var viaAsset = USE_ASSETS && ASSET_CODES.indexOf(u.code) >= 0;

  if (!viaAsset) {
    // One pass. reduceRegions handles empty zones, returning zeros.
    reduceToDrive(zonesOf(u), u);
    return;
  }

  var parts = splitZones(zonesOf(u));
  if (STAGE === 1) {
    Export.table.toAsset({
      collection:  parts.full,
      description: 'zones_' + RUN_TAG + '_' + slug(u.code) +
                   (u.cc ? '_' + u.cc : ''),
      assetId:     assetOf(u)
    });
    // The empty-zone rows never reach the asset, so their zeros are written
    // here, in the pass that already knows which rows they are.
    var zn = nameOf(u) + 'zero';
    Export.table.toDrive(parts.empty.map(function(f){
      return f.setGeometry(null) }), zn, 'TBA2', zn);
  } else {
    reduceToDrive(ee.FeatureCollection(assetOf(u)), u);
  }
});

print(!USE_ASSETS
  ? 'One pass, straight to Drive. CSVs -> 4_typology/1_data/geeOut/.'
  : (STAGE === 1
     ? 'STAGE 1 for ' + ASSET_CODES + ': run, wait, then set STAGE = 2.'
     : 'STAGE 2: CSVs -> 4_typology/1_data/geeOut/.'));

/**** ---------------------------------------------------------------------
IF A TASK FAILS

It will be an S013 country, and only that one — the five European aquifers
here are 16,000-41,000 km² and are not at risk. Options, in order:

  1. USE_ASSETS = true (ASSET_CODES already lists S013): run STAGE 1, wait,
     then STAGE 2. This is what breaks the per-tile recomputation of the
     zone geometry, and it is the reason the option exists at all.
  2. Add S021 / S022 to SIDE_SPLIT_CODES if they turn out to need it — they
     are the next largest at 1.2M and 0.44M km².
  3. For the offending country only, raise SIMPLIFY_M to 1000. The zone edge
     then moves by up to 1 km against a 10 km buffer; note it if you use
     that row.
  4. Last resort: drop the affected rows and report 10 km coverage as
     590/594 rather than forcing it. Missing rows are visible — they simply
     keep their segment-wide values, which is the published behaviour.

VERIFY WHEN DONE
  Rscript 4_typology/gee/missing_dyad_codes.R
should report 10k complete. It diffs keys, not filenames, so it will also
confirm the 4 EU095 - EU098 duplicates are consistent rather than conflicting.
---------------------------------------------------------------------- ****/
