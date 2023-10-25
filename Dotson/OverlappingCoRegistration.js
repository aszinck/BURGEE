/*
This script performs the overlap co-registration of the already CryoSat-2 co-registered REMA strips.
The overlap co-registration is done with respect to a median DEM.

How to run the script:
1. Copy-paste to the code editor on Google Earth Engine
2. Optional: Change data paths in line 35 and 37 in case you want to use your own co-registered DEMs
3. Click run

*/



////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////


// Define period of interest (has to be the median DEM period for which the strips are co-registered to)
var startEpoch = '2010-07-01';
var endEpoch = '2011-07-01';

var startStrip = '2010-07-01';
var endStrip = '2011-07-01';

// Save corrected strips (or not)
var saveToAsset = true;

var overlapCriteria = 25; //%

var allowedTimeSpan = 1;

var outlierCriteria = 150; //meters

var REMA = ee.ImageCollection("users/azcryosphere/DotsonIceShelf/CorrectedREMA_lat60_1month").filterDate(startStrip,endStrip).select(['elevation']);

var REMAint = ee.ImageCollection("users/azcryosphere/DotsonIceShelf/CorrectedREMA_lat60_1month").filterDate(startEpoch,endEpoch).select(['elevation']);


////////////////////////////////////////////////////////////////////
///////////////////////// LOAD DATA ETC. ///////////////////////////
////////////////////////////////////////////////////////////////////

// Initialize map
var AOI =
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-113.69774603371353, -74.15026449134105],
          [-114.05588615686756, -74.24208712637689],
          [-114.04998838474896, -74.35850363192631],
          [-113.95043158058853, -74.47699563098092],
          [-114.3197178806576, -74.72446559450822],
          [-113.94916075340197, -74.96678352337307],
          [-111.93241655226473, -74.89931675089213],
          [-111.47635568874242, -74.83965919930486],
          [-111.53343939308853, -74.75387836080724],
          [-111.35925522843885, -74.45916264740427],
          [-111.58837103371353, -74.34709885255396],
          [-111.74217962746353, -74.22211815622677]]])

Map.setOptions('SATELLITE');
Map.centerObject(AOI,7.5);


// Load ice shelf vector
var IS = ee.FeatureCollection('users/earthmapps/IceShelves').filterMetadata('SURFACE','equals','ice shelf').filterBounds(AOI);
Map.addLayer(IS,{color:'red'},'ice shelf polygon',false);


// Require external functions
var sui  = require('users/earthmapps/public:ui.js');
var palettes = require('users/gena/packages:palettes');
var batch = require('users/fitoprincipe/geetools:batch');

// Use the following line or change to your own directory using the coRegistrationTools script in the modules/ repository
var funcs = require('users/azcryosphere/functions:modules/coRegistrationTools');


// Define colormaps
var elevVis = {min:0, max:100, palette:'450256,3B1C8C,21908D,5AC865,F9E721,white'};


if(REMA.size().getInfo() == 0){
  print('Empty collection');
} else {
var countREMA = REMAint.count().rename('count');
var countMask = countREMA.mask().updateMask(countREMA.select('count').gte(2));

Map.addLayer(countMask,{},'countMask');

var overlapREMA = REMAint.median().updateMask(countMask);

Map.addLayer(overlapREMA,elevVis,'Overlap');


REMA = REMA.map(funcs.percentOverlap(overlapREMA));

var REMAfulfill = REMA.filterMetadata('OverlapPercent','not_less_than',overlapCriteria);
print('REMA fulfill',REMAfulfill);
var REMAnofulfill = REMA.filterMetadata('OverlapPercent','less_than',overlapCriteria);
print('REMA do not fulfill',REMAnofulfill);


var corrREMA = REMAfulfill.map(funcs.tiltBiasCorrectionOverlap(overlapREMA));




corrREMA = corrREMA.map(funcs.finalCorrectionOverlap);
corrREMA = corrREMA.merge(REMAnofulfill);
print('Corrected REMA',corrREMA);
Map.addLayer(REMA,elevVis,'REMA');
Map.addLayer(corrREMA,elevVis,'Corrected REMA');


var exportAOI =
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-114.64075515459402, -74.07451930170846],
          [-114.64075515459402, -75.19479398799669],
          [-110.90540359209402, -75.19479398799669],
          [-110.90540359209402, -74.07451930170846]]], null, false)


if (saveToAsset === true) {
batch.Download.ImageCollection.toAsset(corrREMA,
  'FullyCorrectedREMA',
  {name:'FullyCorrectedREMA',
  scale:10,
  region: exportAOI
})}



Map.addLayer(correctedREMA,elevVis,'Corrected REMA');



var colorleg = ui.Panel({
    style: {
        width: '300px',
        position: 'bottom-right',
        padding: '8px 15px'
    }
});

colorleg.widgets().set(0, sui.colorbar(elevVis,'elevation [m]'));
Map.add(colorleg);
}
