/*
This script calculates the eulerian elevation change

How to run the script:
1. Copy-paste to the code editor on Google Earth Engine
2. Optional: Change data path in line 75 in case you want to use your own co-registered DEMs
3. Click run

Line 126-128 contains more changable parameters controlling the period in which one wants to calculate the Eulerian elevation trend
*/


////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////


// Define period of interest
var startEpoch = '2010-07-01';
var endEpoch = '2018-07-01';

// Elevation change outlier criteria
var maxdhdt = 15; // m/yr

// Max elevation difference compared to the REMA Mosaic
var maxElevDiff = 20; // m

// Save as GEE asset or not
var saveToAsset = true;


////////////////////////////////////////////////////////////////////
///////////////////////// LOAD DATA ETC. ///////////////////////////
////////////////////////////////////////////////////////////////////

// Initialize map

var AOI =
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-113.44910870512868, -74.11557313306545],
          [-114.20716534575368, -74.22645088586815],
          [-114.42689190825368, -74.56335162225412],
          [-114.83338604887868, -74.7638440294105],
          [-115.29481183012868, -75.10931019334664],
          [-114.64661847075368, -75.30844553634742],
          [-113.16346417387868, -75.15441480985159],
          [-111.44959698637868, -74.99596189750002],
          [-111.21888409575368, -74.80421232916828],
          [-111.03211651762868, -74.41941216140597],
          [-111.71326886137868, -74.16661095276221]]]);

Map.setOptions('SATELLITE');
Map.centerObject(AOI,7.5);

// Load ice shelf vector
var IS = ee.FeatureCollection('users/earthmapps/IceShelves').filterMetadata('SURFACE','equals','ice shelf').filterBounds(AOI);
Map.addLayer(IS,{color:'red'},'ice shelf polygon',false);


// Load color maps etc.
var sui  = require('users/earthmapps/public:ui.js');
var palettes = require('users/gena/packages:palettes');

// Use the following line or change to your own directory using the coRegistrationTools script in the modules/ repository
var funcs = require('users/azcryosphere/functions:modules/eulerianTools');

var elevVis = {min:0, max:150, palette:'450256,3B1C8C,21908D,5AC865,F9E721,white'};

var vmin = 10;
var trendVis = {min: -vmin, max: vmin, palette:palettes.cmocean.Balance[7].reverse()};

// Use this imageCollection of fully correted REMA strips or use your own DEMs
var REMA = ee.ImageCollection("users/azcryosphere/DotsonIceShelf/REMADotsonFullCorrection_lat60_1month").filterBounds(AOI).filterDate(startEpoch,endEpoch).select(['elevation']);


////////////////////////////////////////////////////////////////////
////////////////////////// DEM FILTERING ///////////////////////////
////////////////////////////////////////////////////////////////////


var medianREMA = REMA.median();

REMA = REMA.map(funcs.filterByMedianREMA(medianREMA,maxElevDiff));

print(REMA);
////////////////////////////////////////////////////////////////////
/////////////////////// CREATE YEARLY DEMS /////////////////////////
////////////////////////////////////////////////////////////////////

REMA = REMA.map(funcs.addTimeBand);


var calculateEulerian = function(startYear,dt){
  var endYear = startYear + dt;
  var emptyIC = ee.ImageCollection([]);
  for (var i= startYear; i <= endYear; i = i +1) {
    var REMAold = REMA.filterDate(i+'-08-01',(i+1)+'-06-01');

    if(REMAold.size().getInfo() == 0){
      print('empty collection');
      var eulerian = ee.Image.constant(0);
    } else {

    REMAold = REMAold.qualityMosaic('elevation');
    REMAold = REMAold.set('MosaicYear',i+'/'+(i+1));

    emptyIC = emptyIC.merge(REMAold);

  }}
  var eulerian = emptyIC.select(['system:time_start','elevation']).reduce(ee.Reducer.linearFit());
  eulerian = funcs.outlierMask(eulerian.select(['scale']),maxdhdt);
  Map.addLayer(eulerian.clip(IS),trendVis,'eulerian ' + startYear + '/' + (startYear+1) + ' - ' + endYear + '/' + (endYear+1));

  return eulerian.set('dt',dt,'startYear',startYear + '/' + (startYear+1), 'endYear', endYear + '/' + (endYear+1));
};


var yearlyEulerian = ee.ImageCollection([]);

////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////

var startEpoch = 2010; // First year to be considered in the trend period (min is 2010 if using the provided imageCollection which means 2010/11)
var endEpoch = 2017; // Last year to be considered in the trend period (the last data in the provided imageCollection is from 2017/18 --> endEpoch 2017)
var dt = 7; // The number of years to calculate the trend over. Max dt is endEpoch-startEpoch

////////////////////////////////////////////////////////////////////

for (var j = startEpoch; j<= endEpoch-dt; j = j+1) {
  var output = calculateEulerian(j,dt);
  yearlyEulerian = yearlyEulerian.merge(ee.ImageCollection.fromImages([output]));
}

print(yearlyEulerian);


var eulerianFullTrend = yearlyEulerian.first().rename('EulerianElevationChange').clip(IS);
Map.addLayer(eulerianFullTrend,trendVis,'Full eulerian elevation trend');



var colorleg = ui.Panel({
    style: {
        width: '300px',
        position: 'bottom-right',
        padding: '8px 15px'
    }
});
// Add legend colorbars to panel
colorleg.widgets().set(0, sui.colorbar(trendVis,'elevation trend [m/yr]'));
Map.add(colorleg);
