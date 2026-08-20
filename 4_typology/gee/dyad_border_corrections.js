/**** ---------------------------------------------------------------------
Dyad-specific border zones — CORRECTION TABLES
Partitioned by aquifer; supports re-running a subset.
---------------------------------------------------------------------------
Standalone and additive. Does not re-export outCntr_Typol / outTBA_Typol /
outB*k_Typol and does not change how dyads are formed in 52_core.R.

Exports the two zones the preferred specification uses:
    dyadB10k_*.csv   (_B1, near 10 km)
    dyadB100k_*.csv  (_B2, far 100 km)
Drop the parts straight into 4_typology/1_data/geeOut/. apply_dyad_zones()
row-binds every matching part and rejects duplicate keys.

TO RE-RUN ONLY WHAT FAILED
--------------------------
1. Download the parts that DID finish into 1_data/geeOut/.
2. Rscript 4_typology/gee/missing_dyad_codes.R
   It compares what is on disk against the full expected key set and prints
   a ready-to-paste ONLY_CODES array.
3. Paste it below, set RUN_TAG to something new (e.g. 'r1'), run.
   RUN_TAG keeps the new files from colliding with the ones already
   downloaded: dyadB10k_r1_00.csv rather than dyadB10k_00.csv. Both match the
   patterns in 51_config.R.

Do NOT try to re-run "chunk 05" as such. The previous version cut chunks by
list position; this one cuts them by aquifer, so the numbering does not
correspond. Re-run by CODE, which is stable.

WHY THE LAST CHUNK HUNG
-----------------------
The previous version chunked with sides.toList(CHUNK, k * CHUNK). On a
computed collection — and a join result is computed — toList with an offset
is O(offset): to return rows 500-593 the server rebuilds the join and walks
past the first 500. Cost grows with chunk index, so the LAST chunk is always
the worst. Chunks 00-04 land, 05 runs for an hour.

Chunks are now cut by AQUIFER CODE. Dyads only pair countries within the same
aquifer, so partitioning aquifers partitions dyads exactly — nothing is split
or lost. Each task filters the segments to its own codes and joins only that
subset: no offset, no traversal, cost proportional to the group alone.
Groups are balanced on predicted row count n(n-1), longest-first.
---------------------------------------------------------------------- ****/

/**** ---------------------------------------------------------------------
Configuration
---------------------------------------------------------------------- ****/

// Empty = every multi-riparian aquifer. Otherwise only these codes, e.g.
//   var ONLY_CODES = ['AF002','AS126','EU024'];
// missing_dyad_codes.R prints this line for you.
var ONLY_CODES = [];

// Appended to the output names so a re-run cannot overwrite or duplicate
// files already downloaded. '' -> dyadB10k_00.csv; 'r1' -> dyadB10k_r1_00.csv
var RUN_TAG = '';

// Export tasks per width. With ONLY_CODES set, drop this to 1 or 2.
var NGROUPS = 8;

// false: zones straight to Drive.
// true : two passes — STAGE 1 writes geometries to assets, STAGE 2 reduces
//        them. Use only if single-stage tasks stall at every NGROUPS;
//        reduceRegions over COMPUTED geometry re-walks the buffer graph for
//        every raster tile, and materialising breaks that.
var USE_ASSETS   = false;
var STAGE        = 1;
var ASSET_DIR    = 'projects/globaltbav2/assets/VectorBuffers';
var ASSET_PREFIX = 'dyadZones_';

// Only the preferred specification's zones. 5 km and 200 km feed the near5
// and far200 variants and are not run; 51_config.R lists only _B1 and _B2.
var WIDTHS = [[10000,'10k'], [100000,'100k']];

var P_CODE = 'code';
var P_AQ   = 'aq_id';
var P_ISO  = 'CC';

// Exactly the columns 52_core.R reads.
var VARS = ['Area','CR3','GW','GW3','IR','IR3','UR'];

var tiles = 4;
var MIN_RIPARIANS = 3;      // 594 dyad-side rows in total

var RUN_ADJACENCY = false;  // already exported successfully

/**** ---------------------------------------------------------------------
Base rasters — same construction as the production script
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
Riparian counts and the aquifer partition

The histogram is pulled client-side once (426 short keys). It supplies both
the multi-riparian filter and the predicted row count per aquifer, n(n-1),
which is what the groups are balanced on.
---------------------------------------------------------------------- ****/
var nRipDict = ee.Dictionary(TBA.aggregate_histogram(P_CODE)).getInfo();

var multi = [];
for (var c in nRipDict) {
  if (nRipDict[c] < MIN_RIPARIANS) continue;
  if (ONLY_CODES.length && ONLY_CODES.indexOf(c) < 0) continue;
  multi.push({code: c, n: nRipDict[c]});
}
if (ONLY_CODES.length) {
  var found = multi.map(function(a){ return a.code });
  var absent = ONLY_CODES.filter(function(c){ return found.indexOf(c) < 0 });
  if (absent.length)
    print('WARNING: ONLY_CODES entries with no multi-riparian segments:',
          absent);
}

// longest-first: the big aquifers are placed while the groups are still empty
multi.sort(function(a, b){ return b.n - a.n; });

var groups = [], loads = [];
for (var g = 0; g < NGROUPS; g++) { groups.push([]); loads.push(0); }
multi.forEach(function(a){
  var rows = a.n * (a.n - 1);
  var m = 0;
  for (var j = 1; j < NGROUPS; j++) if (loads[j] < loads[m]) m = j;
  groups[m].push(a.code);
  loads[m] += rows;
});

var totalRows = 0;
loads.forEach(function(x){ totalRows += x; });
print('aquifers in this run', multi.length,
      ' dyad-side rows', totalRows,
      ' rows per group', loads);
if (!multi.length)
  print('NOTHING TO RUN: ONLY_CODES matched no multi-riparian aquifer.');

var withN = TBA.map(function(f){
  return f.set('__n', ee.Dictionary(nRipDict).getNumber(f.get(P_CODE)))
})

var pairFilter = ee.Filter.and(
  ee.Filter.equals({leftField: P_CODE,  rightField: P_CODE}),
  ee.Filter.notEquals({leftField: P_ISO, rightField: P_ISO})
)

// Join WITHIN this group's aquifers only.
function sidesOfGroup(codes) {
  var sub = withN.filter(ee.Filter.inList(P_CODE, codes));
  return ee.Join.inner('a','b').apply(sub, sub, pairFilter);
}

function pad(k){ return ('0' + k).slice(-2) }
function partName(tag, k){
  return 'dyadB' + tag + (RUN_TAG ? '_' + RUN_TAG : '') + '_' + pad(k);
}
function assetOf(tag, k){
  return ASSET_DIR + '/' + ASSET_PREFIX + tag +
         (RUN_TAG ? '_' + RUN_TAG : '') + '_' + pad(k);
}

function keysOf(A, B){
  return {
    code:     A.get(P_CODE),
    aq_id:    A.get(P_AQ),
    CC:       A.get(P_ISO),
    CC_other: B.get(P_ISO),
    n_rip:    A.get('__n')
  }
}

/**** ---------------------------------------------------------------------
Zone geometry

zone = segment_A ∩ segment_B.buffer(d)

A pixel in segment A within d of segment B is within d of the A-B border:
reaching B from A means crossing it, and the nearest point of B from inside A
lies on the shared boundary. Clipping to the aquifer is automatic — both
segments already are the aquifer.

Tolerances scale with width: 10 km -> simplify 500 m, maxError 100 m;
100 km -> 2 km, 500 m. Both below the 1 km reduction grid, and buffer cost
tracks vertex count.
---------------------------------------------------------------------- ****/
function zonesOf(codes, d) {
  var simp  = Math.max(500, d / 50);
  var bufME = Math.max(100, d / 200);
  return sidesOfGroup(codes).map(function(row){
    var A = ee.Feature(row.get('a'));
    var B = ee.Feature(row.get('b'));
    var zone = A.geometry().intersection(
      B.geometry().simplify(simp).buffer(d, bufME), bufME);
    return ee.Feature(zone, keysOf(A, B));
  });
}

function reduceAndExport(fc, tag, k) {
  var out = IM.reduceRegions({
    collection: fc,
    reducer:    ee.Reducer.sum(),
    scale:      1000,
    tileScale:  tiles
  }).map(function(f){ return f.setGeometry(null) });
  var nm = partName(tag, k);
  Export.table.toDrive(out, nm, 'TBA2', nm);
}

/**** ---------------------------------------------------------------------
Export
---------------------------------------------------------------------- ****/
WIDTHS.forEach(function(w){
  var d = w[0], tag = w[1];
  for (var k = 0; k < NGROUPS; k++) {
    if (!groups[k].length) continue;
    if (!USE_ASSETS) {
      reduceAndExport(zonesOf(groups[k], d), tag, k);
    } else if (STAGE === 1) {
      Export.table.toAsset({
        collection:  zonesOf(groups[k], d),
        description: 'zones_' + tag + '_' + pad(k),
        assetId:     assetOf(tag, k)
      });
    } else {
      reduceAndExport(ee.FeatureCollection(assetOf(tag, k)), tag, k);
    }
  }
});

/**** ---------------------------------------------------------------------
Adjacency / QA — already exported successfully
---------------------------------------------------------------------- ****/
if (RUN_ADJACENCY) {
  var adjParts = [];
  for (var k2 = 0; k2 < NGROUPS; k2++) {
    if (!groups[k2].length) continue;
    adjParts.push(sidesOfGroup(groups[k2]).map(function(row){
      var A = ee.Feature(row.get('a'));
      var B = ee.Feature(row.get('b'));
      return ee.Feature(null, keysOf(A, B))
        .set('gap_m', A.geometry().distance(B.geometry(), 100));
    }));
  }
  var an = 'dyadAdjacency' + (RUN_TAG ? '_' + RUN_TAG : '');
  Export.table.toDrive(ee.FeatureCollection(adjParts).flatten(),
                       an, 'TBA2', an);
}

/**** ---------------------------------------------------------------------
CHECKS once everything is in geeOut/

  1. Rscript gee/missing_dyad_codes.R reports nothing missing. It counts the
     expected 594 (code, CC, CC_other) keys per width from outTBA_Typol.csv
     and diffs them against what is on disk — so it catches a part that was
     never downloaded, not only one that never ran.
  2. No duplicate key. apply_dyad_zones() stops rather than summing, which is
     the failure mode a partial re-run invites.
  3. Every value <= the matching outB{10k,100k}_Typol value for the same
     (code, CC): the dyad border is a subset of all land borders.
---------------------------------------------------------------------- ****/
