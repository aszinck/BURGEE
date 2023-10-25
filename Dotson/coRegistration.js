/*
This script performs the initial co-registration of the REMA 2m strips over the Dotson Ice Shelf with respect to CryoSat-2

How to run the script:
Copy-paste to the code editor on Google Earth Engine and click run
*/


////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////


// Define period of interest
var startEpoch = '2010-07-01';
var endEpoch = '2011-07-01';

// Save corrected strips (or not)
var saveToAsset = true;

// Define period allowed for co-registration
var allowedTimeSpan = 1; //months

// Define the max difference allowed between a CS2 and REMA elevation
var outlierCriteria = 100; // meters

// Define the number of matchups required to perform the co-registation
var numMatchups = 80; //points

// Define latitude criteria for CryoSat distribution
var latCriteria = 60; //km

// Define longitude criteria for CryoSat distribution
var lonCriteria = 10; //km

// Max difference from REMA Mosaic criteria
var elevMosaic = 30; //m

// Mean sea level pressure
var SLP = 1013;



////////////////////////////////////////////////////////////////////
///////////////////////// LOAD DATA ETC. ///////////////////////////
////////////////////////////////////////////////////////////////////

// Initilize map
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


// Load ice shelf geometry
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

// Import REMA 2-m strips
var REMA = ee.ImageCollection("UMN/PGC/REMA/V1/2m")

// Filter REMA based on ice shelf and period of interest
REMA = REMA.filterBounds(IS).filterDate(startEpoch,endEpoch).select(['elevation']);

// Exclude strips from GeoEye-1 due to stripes
REMA = REMA.filterMetadata('platform1','not_equals','GE01');

print('Initial REMA imageCollection',REMA);

print('Initial number of REMA strips:',REMA.size());

// Import CryoSat-2 data
var CS2 = ee.FeatureCollection("users/azcryosphere/CryoSat/CryoSat2_FullDotson");

// Add proper timestamp to the CS-2 data and filter based on period of interest
CS2 = CS2.map(funcs.CS2systemtimestart);
CS2 = CS2.filterDate(ee.Date(startEpoch).advance(-allowedTimeSpan, 'month'),ee.Date(endEpoch).advance(allowedTimeSpan, 'month'));


////////////////////////////////////////////////////////////////////
////////////////////// REMA MOSAIC FILTERING ///////////////////////
////////////////////////////////////////////////////////////////////

// Import REMA mosaic
var REMAMosaic = ee.Image("UMN/PGC/REMA/V1_1/8m")

// Mask out points on REMA strips that differ more than *elevMosaic* from the REMA Mosaic
REMA = REMA.map(funcs.filterByREMAMosaic(elevMosaic,REMAMosaic));



////////////////////////////////////////////////////////////////////
/////////////////////// DYNAMIC CORRECTIONS ////////////////////////
////////////////////////////////////////////////////////////////////

// Import geoid data
var EGM = ee.Image("users/azcryosphere/GlaciologicalData/egm2008-1")

// Perform all dynamic corrections on CS-2
CS2 = CS2.map(funcs.CryoSatDynamicCorrections(IS,SLP,EGM,ee.Date(startEpoch).advance(-allowedTimeSpan, 'month'),ee.Date(endEpoch).advance(allowedTimeSpan, 'month')));
CS2 = CS2.filterMetadata('NCEP','equals','present');
CS2 = CS2.select(['h','system:time_start']);

// Perform all dynamic corrections on REMA strips
REMA = REMA.map(funcs.REMAdynamicCorrectionsV2(IS,SLP,EGM));



////////////////////////////////////////////////////////////////////
////////////////////// MASKING AND FILTERING ///////////////////////
////////////////////////////////////////////////////////////////////



// Filter by number of CS-2 points on one strip
REMA = REMA.map(funcs.filterByCryoSat(CS2,allowedTimeSpan,outlierCriteria));
REMA = REMA.filterMetadata('matchups','greater_than',numMatchups);


// Filter by CS-2 spatial distribution
REMA = REMA.map(funcs.filterByCryoSatDistribution(CS2,allowedTimeSpan));
REMA = REMA.filterMetadata('latDiff','greater_than',latCriteria)
          .filterMetadata('lonDiff','greater_than',lonCriteria);



// ////////////////////////////////////////////////////////////////////
// ///////////////////////// CO-REGISTRATION //////////////////////////
// ////////////////////////////////////////////////////////////////////

// Perform co-registration
REMA = REMA.map(funcs.tiltBiasCorrection(CS2,allowedTimeSpan,outlierCriteria));

// Apply co-registration
var correctedREMA = REMA.map(funcs.finalCorrection);
print(correctedREMA);

print('Final number of REMA strips:',correctedREMA.size());

// Define the export region
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
batch.Download.ImageCollection.toAsset(correctedREMA,
  'CorrectedREMA',
  {name:'CorrectedREMA',
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
