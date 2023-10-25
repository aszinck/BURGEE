////////////////////////////////////////////////////////////////////
////////////////////////// GENERAL IMPORTS /////////////////////////
////////////////////////////////////////////////////////////////////

// Import geoid 
var EGM = ee.Image("users/azcryosphere/GlaciologicalData/egm2008-1");

// Import the 2-meter REMA strips
var REMA2 = ee.ImageCollection("UMN/PGC/REMA/V1/2m");

// Import the 8-meter REMA strips
var REMA8 = ee.ImageCollection("UMN/PGC/REMA/V1/8m");

// Define an export Area Of Interest which should match *AOI* and *GL*
// *exportAOI* should not be bigger than *AOI* or *GL*
// With the *exportAOI* there is the option to exclude ice shelves that are only partially covered by *AOI* or *GL*
var exportAOI = 
    /* color: #ffc82d */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[80.15439569184237, -65.36484608696702],
          [80.14803107096165, -66.19769756691068],
          [80.14451285865252, -67.45202965417329],
          [80.108959300282, -68.99988655223264],
          [82.63473278088091, -68.98914567825493],
          [86.58840262421239, -68.97919532784852],
          [89.96884408460464, -68.98561284699505],
          [92.14275242965711, -68.97919690909012],
          [95.65396692933058, -68.98881719937873],
          [99.25689649814488, -68.98945799929567],
          [101.4359596078105, -68.97855471948219],
          [104.93026034628097, -68.99053531690836],
          [104.92936131658143, -66.6780588552092],
          [104.95132916455356, -66.06945092077808],
          [104.96258173559015, -65.21241529751124],
          [104.99830314617529, -64.43746019500979],
          [101.02201997349755, -64.75543054728566],
          [97.78031249610356, -64.91423366278936],
          [95.47874059404481, -65.00474989881944],
          [93.52777401333081, -65.04598195708493],
          [91.21604537587478, -65.05502034206171],
          [89.71851432089734, -65.01690512215798],
          [88.39446935250751, -64.97539507866908],
          [86.38872139508106, -64.90615432518399],
          [83.41771521628479, -64.73731463779932],
          [80.1545019532602, -64.44921904029368]]]);

// Load the ocean mask
var oceanmask = ee.Image("users/aszinck/GlaciologicalData/BMoceanmask");


////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////

// Define region of interest
var region = 'ZoneWilkes1';           // Possible regions: ZoneAmundsen1, ZoneAmundsen2, ZoneAmundsen3, ZoneBellinghausen1, ZOneBellinghausen2, ZoneVictoria1, ZoneVictoria2, ZoneWilkes1, ZoneWilkes2, ZoneLarsen 

// Define period of interest
// OBS! The possible period length possible to handle in one run depends on the number of REMA strips within that period. Shorter periods may be necessary sometimes.
var startEpoch = '2010-07-01'; 
var endEpoch = '2011-07-01'; 

// Save corrected strips (or not)
var saveToAsset = true;

// Define period allowed for co-registration 
var allowedTimeSpan = 1; // [months]           The timedifference allowed between a REMA strips and the CryoSat-2 data

// Define the number of matchups and coverage required to perform the co-registation
var numMatchups = 80; // [points]
var numPercentageCS2coverage = 0.05; // [%]
var numMatchupsLow = 5; // [points]

// Define latitude criteria for CryoSat distribution
var latCriteria = 0.75; // [%]

// Define longitude criteria for CryoSat distribution
var lonCriteria = 0.75; // [%]

// Max difference from REMA Mosaic criteria
var elevMosaic = 30; // [m] 

// Mean sea level pressure
var SLP = 1013; // [hPa]



////////////////////////////////////////////////////////////////////
///////////////////////// LOAD DATA ETC. ///////////////////////////
////////////////////////////////////////////////////////////////////

// Load AOI based on region
var AOI = ee.FeatureCollection('users/aszinck/GlaciologicalData/regionsCollection').filterMetadata('region','equals',region).geometry();
Map.addLayer(AOI,{color:'black'},'AOI',false);

Map.setOptions('SATELLITE');
Map.centerObject(AOI,6);


print('This script performs the initial co-registration \nof the ice shelves in the '+region+' region.');

print('Make sure to change exportAOI to cover the right \nice shelves within the region.');


// Load ice shelf vector
var IS = ee.FeatureCollection('users/azcryosphere/GlaciologicalData/FilteredIS').filterBounds(exportAOI);
Map.addLayer(IS,{color:'red'},'ice shelf polygon',false);

// Load grounding line
var GL = ee.FeatureCollection('users/aszinck/GlaciologicalData/'+region+'Region').filterMetadata('IS','equals',0).select('IS').filterBounds(AOI);
Map.addLayer(GL,{color:'green'},'GL polygon',false);

// Require external functions
var sui  = require('users/earthmapps/public:ui.js'); // Colorbar tools
var palettes = require('users/gena/packages:palettes'); // Different colormaps
var tools = require('users/fitoprincipe/geetools:tools');  // Different batch export tools
var funcs = require('users/azcryosphere/functions:modules/coRegistrationUniversalTools'); // co-registration tools



// Define colormaps
var elevVis = {min:0, max:100, palette:'450256,3B1C8C,21908D,5AC865,F9E721,white'};

// Merge the 2- and 8-meter strips into one imageCollection
var REMA = REMA2.merge(REMA8);

// Filter REMA based on ice shelf and period of interest
REMA = REMA.filterBounds(IS).filterDate(startEpoch,endEpoch).select(['elevation']);

// Exclude strips from GeoEye-1 due to stripes
REMA = REMA.filterMetadata('platform1','not_equals','GE01');

print('REMA',REMA);

Map.addLayer(REMA,elevVis,'Initial strips');

print('Initial number of REMA strips:',REMA.size());



// Add proper timestamp to the CryoSat-2 data and filter based on period of interest
var CS2 = ee.FeatureCollection('users/aszinck/CryoSat/CryoSat2_'+region).filterBounds(REMA.geometry());
CS2 = CS2.map(funcs.CS2systemtimestart);
CS2 = CS2.filterDate(ee.Date(startEpoch).advance(-allowedTimeSpan, 'month'),ee.Date(endEpoch).advance(allowedTimeSpan, 'month'));

// Map.addLayer(CS2,{},'CS2')

////////////////////////////////////////////////////////////////////
////////////////////// REMA MOSAIC FILTERING ///////////////////////
////////////////////////////////////////////////////////////////////


// Mask out points on REMA strips that differ more than *elevMosaic* from the REMA Mosaic
var REMAMosaic = ee.Image('users/aszinck/REMAmosaics/'+region+'REMAdemMosaic').updateMask(oceanmask.mask());
REMA = REMA.map(funcs.filterByREMAMosaic(elevMosaic,REMAMosaic));

// print('filterByREMAMosaic',REMA)



////////////////////////////////////////////////////////////////////
////////////////// DYNAMIC AND STATIC CORRECTIONS //////////////////
////////////////////////////////////////////////////////////////////

/*
Dynamic and static corrections: 
- Tides
- Inverse barometer effect
- Mean dynamic topography
- Geoid
*/



// Perform all dynamic corrections on CS-2
CS2 = CS2.map(funcs.CryoSatDynamicCorrectionsAlpha(SLP,EGM,ee.Date(startEpoch).advance(-allowedTimeSpan, 'month'),ee.Date(endEpoch).advance(allowedTimeSpan, 'month'),region));
CS2 = CS2.filterMetadata('NCEP','equals','present');
CS2 = CS2.select(['h','system:time_start']);

// print(CS2.first())

  

// Perform all corrections on REMA strips
REMA = REMA.map(funcs.REMAdynamicCorrectionsAlpha(GL,SLP,EGM,region));


////////////////////////////////////////////////////////////////////
////////////////////// MASKING AND FILTERING ///////////////////////
////////////////////////////////////////////////////////////////////


// Filter by percentage coverage of CS-2 points at a 500 m res on one strip
REMA = REMA.map(funcs.filterByCryoSatPercentage(CS2,allowedTimeSpan)); // Adds information how high the CS-2 coverage is in percentage with respect to the REMA strips area
REMA = REMA.map(funcs.filterByCryoSat(CS2,allowedTimeSpan)); // Adds information about how many CS-2 points there are available within one REMA strips

REMA = REMA.filterMetadata('matchups','not_less_than',numMatchupsLow); // Filter out all strips with less than *numMatchupsLow* CS-2 points

// The following steps ensures that either *numPercentageCS2coverage* or *numMatchups* is fulfilled
var REMAp = REMA.filterMetadata('percentageCS2coverage','not_less_than',numPercentageCS2coverage);
var REMAnp = REMA.filterMetadata('percentageCS2coverage','less_than',numPercentageCS2coverage);

print('REMA percentage fulfill',REMAp);
REMAnp = REMAnp.filterMetadata('matchups','not_less_than',numMatchups);
print('REMA matchup only fulfill',REMAnp);

var REMA = REMAp.merge(REMAnp);
print('After CS2 percentage and matchups',REMA);



// Filter by CS-2 spatial distribution
REMA = REMA.map(funcs.filterByCryoSatCoverage(CS2,allowedTimeSpan));
REMA = REMA.filterMetadata('latFrac','greater_than',latCriteria)
          .filterMetadata('lonFrac','greater_than',lonCriteria); 

print('filterByCryoSatDistribution',REMA);



// ////////////////////////////////////////////////////////////////////
// ///////////////////////// CO-REGISTRATION //////////////////////////
// ////////////////////////////////////////////////////////////////////

// Perform co-registration
REMA = REMA.map(funcs.tiltBiasCorrection(CS2,allowedTimeSpan));

// Apply co-registration
var correctedREMA = REMA.map(funcs.finalCorrection);
// print(correctedREMA);




////////////////////////////////////////////////////////////////////
//////////////////////////// IS MASKING ////////////////////////////
////////////////////////////////////////////////////////////////////

// Discard strips that do not have data on any of the ice shelves of interest

function clipCollection(image){
  var mask = ee.Image.constant(1).clip(IS).mask();
  image = image.updateMask(mask);
  var areaImage = ee.Image.constant(1).updateMask(image.mask());
  var area = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: IS.geometry(),
    scale: 500
    });

  return image.set('Area',ee.Number(area.get('constant')));
}



correctedREMA = correctedREMA.map(clipCollection);

correctedREMA = correctedREMA.filter(ee.Filter.metadata('Area','greater_than',1));

print('After IS masking',correctedREMA);




print('Final number of REMA strips:',correctedREMA.size());


var downloadAll = function(collection, assetFolder, options) {
  // var root = ee.data.getAssetRoots()[0]['id']
  // var root = 'aszinck'
  // var folder = assetFolder
  // var assetFolder = root+'/'+folder+'/'
  // if (folder !== null && folder !== undefined) {
  //   var assetFolder = root+'/'+folder+'/'
  // } else {
  //   var assetFolder = root+'/'
  // }
  
  var defaults = {
      name: null,
      scale: 1000,
      maxPixels: 1e13,
      region: null
    };
    
  var opt = tools.get_options(defaults, options);
  var n = collection.size().getInfo();
    
  var colList = collection.toList(n);
  
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = img.id().getInfo() || 'image_'+i.toString();
    // var region = opt.region || img.geometry(1).bounds(1).getInfo()["coordinates"];
    var region = opt.region || img.geometry(1);
    var assetId = assetFolder+'/'+id;
    
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: assetId,
      region: region,
      scale: opt.scale,
      maxPixels: opt.maxPixels});
  }
};



// Save the REMA strips to the inagecollection CorrectedREMA within the *region* folder
if (saveToAsset === true) { 
downloadAll(correctedREMA,
  region+'/CorrectedREMA', 
  {name:'CorrectedREMA',
  scale:50
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
