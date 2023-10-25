

////////////////////////////////////////////////////////////////////
/////////////////////// CHANGABLE PARAMETERS ///////////////////////
////////////////////////////////////////////////////////////////////

var region = 'ZoneWilkes1'; 


// Define period of interest (has to be the DEM period)
var startEpoch = '2010-07-01'; 
var endEpoch = '2011-07-01'; 

var startStrip = '2010-07-01'; 
var endStrip = '2011-07-01'; 

// Save corrected strips (or not)
var saveToAsset = true;

var overlapCriteria = 25; //%


var AOI = ee.FeatureCollection('users/aszinck/GlaciologicalData/regionsCollection').filterMetadata('region','equals',region).geometry();
Map.addLayer(AOI,{color:'black'},'AOI',false);

Map.setOptions('SATELLITE');
Map.centerObject(AOI,5);

print('This script performs the second co-registration \nof the the ice shelves in the '+region+' region.');



////////////////////////////////////////////////////////////////////
///////////////////////// LOAD DATA ETC. ///////////////////////////
////////////////////////////////////////////////////////////////////


// Load ice shelf vector
var IS = ee.FeatureCollection('users/azcryosphere/BedMachineGroundingLine_50m').filterMetadata('IS','equals',1).select('IS').filterBounds(AOI);
Map.addLayer(IS,{color:'red'},'ice shelf polygon',false);


// Require external functions
var sui  = require('users/earthmapps/public:ui.js');
var palettes = require('users/gena/packages:palettes');
var funcs = require('users/azcryosphere/functions:modules/coRegistrationUniversalTools');
var batch = require('users/fitoprincipe/geetools:batch');
var tools = require('users/fitoprincipe/geetools:tools'); 


// Define colormaps
var elevVis = {min:0, max:100, palette:'450256,3B1C8C,21908D,5AC865,F9E721,white'};

var REMA = ee.ImageCollection('users/azcryosphere/'+region+'/CorrectedREMA').filterDate(startStrip,endStrip).select(['elevation']);
print('Old corrected',REMA);

var REMAint = ee.ImageCollection('users/azcryosphere/'+region+'/CorrectedREMA').filterDate(startEpoch,endEpoch).select(['elevation']);

Map.addLayer(REMAint,elevVis,'Orig REMA');



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
print('REMA fulfill',REMAfulfill.size());
var REMAnofulfill = REMA.filterMetadata('OverlapPercent','less_than',overlapCriteria);
print('REMA do not fulfill',REMAnofulfill.size());


var corrREMA = REMAfulfill.map(funcs.tiltBiasCorrectionOverlap(overlapREMA));



corrREMA = corrREMA.map(funcs.finalCorrectionOverlap);
corrREMA = corrREMA.merge(REMAnofulfill);
print('Corrected REMA',corrREMA);
Map.addLayer(REMAnofulfill,elevVis,'REMA no fulfill');
Map.addLayer(corrREMA,elevVis,'Corrected REMA');


var downloadAll = function(collection, assetFolder, options) {
  
  var defaults = {
      name: null,
      scale: 1000,
      maxPixels: 1e13,
      region: null
    }
    
  var opt = tools.get_options(defaults, options)
  var n = collection.size().getInfo();
    
  var colList = collection.toList(n);
  
  for (var i = 0; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = img.id().getInfo() || 'image_'+i.toString();
    // var region = opt.region || img.geometry(1).bounds(1).getInfo()["coordinates"];
    var region = opt.region || img.geometry(1);
    var assetId = assetFolder+'/'+id
    
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: assetId,
      region: region,
      scale: opt.scale,
      maxPixels: opt.maxPixels})
  }
}




if (saveToAsset === true) { 
downloadAll(corrREMA,
  region+'/FullyCorrectedREMA', 
  {name:'FullyCorrectedREMA',
  scale:50
})}




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