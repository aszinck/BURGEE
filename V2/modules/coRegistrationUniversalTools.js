/////////////////////////////////////////////////////////////////////////////
///////////////////////// CO-REGISTRATION TOOLS /////////////////////////////
/////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////
////////////////////// ADD system:time_start TO CS2 /////////////////////////
/////////////////////////////////////////////////////////////////////////////
 
 
var CS2systemtimestart = function (feature) {
  
  /* Add a proper system:time_start field to the CS-2 featureCollection
  */
  
  var dateC = ee.Date.parse('YYYY-MM-dd HH:mm:ss',ee.String(feature.get('timestamp')));
  return feature.set('system:time_start', dateC.millis());
};

exports.CS2systemtimestart = CS2systemtimestart;




/////////////////////////////////////////////////////////////////////////////
////////// ATMOSPHERIC AND OCEANIC CORRECTIONS FOR CS2 AND REMA /////////////
/////////////////////////////////////////////////////////////////////////////



var roundDateTo6h = function(date){
  
  /* Round any date to nearest 6h interval
  */
  
  var roundTo = ee.Number(60*60*6);
  var seconds = date.difference(date.update(date.get('year'),date.get('month'),date.get('day'),0,0,0),'seconds');
  var rounding = ((seconds.add(roundTo.divide(ee.Number(2)))).divide(roundTo).floor()).multiply(roundTo);
  return date.advance(rounding.subtract(seconds),'seconds');
};






var CryoSatDynamicCorrectionsAlpha = function(SLP,EGM,startEpoch,endEpoch,region) {
  
  /* Correct a featureCollection for tide, MDT, IBE and geoid
  *
  *  IBE: I use a correction of 0.9948 cm/hPa (http://www.altimetry.info/radar-altimetry-tutorial/data-flow/data-processing/geophysical-corrections/inverse-barometer/)
  *  'NCEP',answer provides information as to whether there is any NCEP data present from the given timestamp or not
  *
  */
  
  return(
  function(feature){

  var alpha = ee.Image('users/aszinck/GlaciologicalData/alpha'+region);
  
  // MDT
  var MDT = ee.Image('users/azcryosphere/GlaciologicalData/DTU15MDT_1dg');
  MDT = MDT.multiply(alpha);
  var MDT_CS = MDT.reduceRegion({reducer:ee.Reducer.mean(), geometry:feature.geometry(), scale:1}).get('b1');
  
  // IBE
  var NCEP = ee.ImageCollection("NCEP_RE/sea_level_pressure").filterDate(feature.get('system:time_start'));
  var sizeNcep = NCEP.size();
  var answer = ee.Algorithms.If(sizeNcep, "present", "NOTpresent");
  var NCEPdate = ee.Algorithms.If(sizeNcep, NCEP.first().divide(ee.Image.constant(100)), ee.Image.constant(0).rename('slp'));
  NCEPdate = ee.Image(NCEPdate).multiply(alpha);
  var pressure = NCEPdate.reduceRegion({reducer:ee.Reducer.mean(), geometry:feature.geometry(), scale:1}).get('slp');
  var IBE_CS = ee.Number(SLP).subtract(pressure);
  IBE_CS = IBE_CS.divide(100).multiply(0.9948);

  // geoid
  var geoid_CS = EGM.reduceRegion({reducer:ee.Reducer.mean(), geometry:feature.geometry()}).get('b1');

  // tide
  var tide = ee.Number(feature.get('h_tide'));
  
  var h = ee.Number(feature.get('h_ell')).subtract(tide).subtract(ee.Number(MDT_CS)).subtract(ee.Number(IBE_CS)).subtract(ee.Number(geoid_CS));
  
  return feature.set('h',h,'MDT',MDT_CS,'IBE',IBE_CS,'geoid',geoid_CS,'NCEP',answer);
})};


exports.CryoSatDynamicCorrectionsAlpha = CryoSatDynamicCorrectionsAlpha;





var REMAdynamicCorrectionsAlpha = function(GL,SLP,EGM,region) {
  
  /* Correct REMA strips for MDT, IBE, tides, and geoid.
  *
  *  IBE: I use a correction of 0.9948 cm/hPa (http://www.altimetry.info/radar-altimetry-tutorial/data-flow/data-processing/geophysical-corrections/inverse-barometer/)
  *
  */
  
  return(
  function(image){
    

  var alpha = ee.Image('users/aszinck/GlaciologicalData/alpha'+region);
      
  var sID = ee.String(image.get('sourceImage1'));
  var mDate = sID.split('_').get(1);
  var newDate = ee.Date.parse('YYYYMMdHHmmss',ee.String(mDate));
  var imageDate = roundDateTo6h(newDate);
  
  // MDT
  var maskREMA = image.mask();
  var MDT = ee.Image('users/azcryosphere/GlaciologicalData/DTU15MDT_1dg');
  MDT = MDT.multiply(alpha);
  
  // IBE
  var NCEP = ee.ImageCollection("NCEP_RE/sea_level_pressure").filterDate('2000-01-01','2022-01-01');
  var meanPressure = ee.Image.constant(SLP).updateMask(maskREMA);
  var NCEPdate = NCEP.filterDate(imageDate).first().divide(ee.Image.constant(100));
  var pressure = NCEPdate.subtract(meanPressure);
  var IBE = pressure.divide(ee.Image.constant(1).clip(GL).mask()).divide(100).multiply(0.9948).multiply(alpha); 
  
  // Tides
  var tide = ee.ImageCollection("users/aszinck/Tides/"+region+"_CATS2008").filterDate(newDate).first();
  tide = tide.divide(ee.Image.constant(1).clip(GL).mask());
  tide = tide.unmask(0,false);
  tide = tide.multiply(alpha);
  
  image = image.subtract(IBE).subtract(MDT).subtract(EGM).subtract(tide).set("system:time_start", image.get("system:time_start"));
  
  return image;
})};

exports.REMAdynamicCorrectionsAlpha = REMAdynamicCorrectionsAlpha;



/////////////////////////////////////////////////////////////////////////////
////////////////////////// MASKING AND FILTERING ////////////////////////////
/////////////////////////////////////////////////////////////////////////////



var filterByCryoSat = function(CS2,allowedTimeSpan){
  
  /* Masking and filtering
  *
  *  1. Mask out points where the abs. elev. diff. between REMA and CS-2 exceeds *outlierCriteria*
  *  2. Count the number of CS-2 points within the REMA strip and add *matchups* to the REMA metadata
  *
  */
  
  return(
  function(image){
  var imageDate = ee.Date(image.get('system:time_start'));
  var startDate = imageDate.advance(-allowedTimeSpan, 'month');
  var endDate = imageDate.advance(allowedTimeSpan, 'month');
  var CS2_image = CS2.filterDate(startDate,endDate).filterBounds(image.geometry());
  CS2_image = CS2_image.reduceToImage(['h'], ee.Reducer.mean()).rename('elevation');
  var maskREMA = image.mask();
  var maskCS = CS2_image.mask();
  var diff = CS2_image.subtract(image).abs();
  diff = diff.updateMask(maskREMA).updateMask(maskCS).rename('diff');
  var reg = image.geometry();
  var scale = 500;
  var constantImage = ee.Image.constant(1).updateMask(maskREMA).updateMask(maskCS).updateMask(diff.select('diff')); //.lte(outlierCriteria));
  var matchups = constantImage.reduceRegion({
    reducer:ee.Reducer.sum(), 
    geometry:reg,
    scale:scale,
    crs:'EPSG:3031',
    tileScale:16})
    .get('constant');

  return image.set({matchups:matchups});
})};

exports.filterByCryoSat = filterByCryoSat;


var filterByCryoSatPercentage = function(CS2,allowedTimeSpan){
  
  /* Masking and filtering
  *
  *  1. Mask out points where the abs. elev. diff. between REMA and CS-2 exceeds *outlierCriteria*
  *  2. Count the number of CS-2 points within the REMA strip and add *matchups* to the REMA metadata
  *
  */
  
  return(
  function(image){
  var imageDate = ee.Date(image.get('system:time_start'));
  var startDate = imageDate.advance(-allowedTimeSpan, 'month');
  var endDate = imageDate.advance(allowedTimeSpan, 'month');
  var CS2_image = CS2.filterDate(startDate,endDate).filterBounds(image.geometry());
  CS2_image = CS2_image.reduceToImage(['h'], ee.Reducer.mean()).rename('elevation');
  var maskREMA = image.mask();
  var maskCS = CS2_image.mask();
  var diff = CS2_image.subtract(image).abs();
  diff = diff.updateMask(maskREMA).updateMask(maskCS).rename('diff');
  var reg = image.geometry();
  var scale = 500;
  var constantImage = ee.Image.constant(1).updateMask(maskREMA).updateMask(maskCS).updateMask(diff.select('diff')); //.lte(outlierCriteria));
  var matchups = constantImage.reduceRegion({
    reducer:ee.Reducer.sum(), 
    geometry:reg,
    scale:scale,
    crs:'EPSG:3031',
    tileScale:16})
    .get('constant');
  var constantREMA = ee.Image.constant(1).updateMask(maskREMA);
  var matchupsREMA = constantREMA.reduceRegion({
    reducer:ee.Reducer.sum(), 
    geometry:reg,
    scale:scale,
    crs:'EPSG:3031',
    tileScale:16})
    .get('constant');
  var percentageCS2coverage = ee.Number(matchups).divide(ee.Number(matchupsREMA));
  image = image.set({matchupsREMA:matchupsREMA});
  image = image.set({matchups:matchups});
  image = image.set({percentageCS2coverage:percentageCS2coverage});
  return image;
})};

exports.filterByCryoSatPercentage = filterByCryoSatPercentage;


var filterByCryoSatDistribution = function(CS2,allowedTimeSpan){
  
  /* Filter by CS-2 spatial distribution
  *
  *  The functions adds to new pieces of metadata:
  *  *latDiff* the distance [km] between the two CS-2 points with the min. and max. latitude
  *  *lonDiff* the distance [km] between the two CS-2 points with the min. and max. longitude
  *
  */
  
  return(
  function(image){
  var imageDate = ee.Date(image.get('system:time_start'));
  var startDate = imageDate.advance(-allowedTimeSpan, 'month');
  var endDate = imageDate.advance(allowedTimeSpan, 'month');
  var CS2_image = CS2.filterDate(startDate,endDate).filterBounds(image.geometry());
  CS2_image = CS2_image.reduceToImage(['h'], ee.Reducer.mean()).rename('elevation');
  var maskREMA = image.mask();
  var maskCS = CS2_image.mask();
  var latitude = ee.Image.pixelLonLat().select('latitude').updateMask(maskREMA).updateMask(maskCS);
  var latMinMax = latitude.reduceRegion({reducer: ee.Reducer.minMax(), geometry: image.geometry(), scale: 500, crs: 'EPSG:3031'});
  var longitude = ee.Image.pixelLonLat().select('longitude').updateMask(maskREMA).updateMask(maskCS);
  var lonMinMax = longitude.reduceRegion({reducer: ee.Reducer.minMax(), geometry: image.geometry(), scale: 500, crs: 'EPSG:3031'});
  var center = image.geometry().centroid({'maxError': 1}).coordinates().get(1);
  var deg2km = ee.Number(40000/360);
  var deg2rad = ee.Number(center).multiply(ee.Number(Math.PI)).divide(ee.Number(180));
  var lon2km = deg2km.multiply(deg2rad.cos());
  var latDiff = (ee.Number(latMinMax.get('latitude_max')).subtract(ee.Number(latMinMax.get('latitude_min')))).multiply(deg2km);
  var lonDiff = (ee.Number(lonMinMax.get('longitude_max')).subtract(ee.Number(lonMinMax.get('longitude_min')))).multiply(lon2km);

  return image.set({latDiff:latDiff, lonDiff:lonDiff});
})};

exports.filterByCryoSatDistribution = filterByCryoSatDistribution;



var filterByCryoSatCoverage = function(CS2,allowedTimeSpan){
  
  /* Filter by CS-2 spatial distribution
  *
  *  The functions adds to new pieces of metadata:
  *  *latDiff* the distance [km] between the two CS-2 points with the min. and max. latitude
  *  *lonDiff* the distance [km] between the two CS-2 points with the min. and max. longitude
  *
  */
  
  return(
  function(image){
  var imageDate = ee.Date(image.get('system:time_start'));
  var startDate = imageDate.advance(-allowedTimeSpan, 'month');
  var endDate = imageDate.advance(allowedTimeSpan, 'month');
  var CS2_image = CS2.filterDate(startDate,endDate).filterBounds(image.geometry());
  CS2_image = CS2_image.reduceToImage(['h'], ee.Reducer.mean()).rename('elevation');
  var maskREMA = image.mask();
  var maskCS = CS2_image.mask();
  
  //// REMA size
  var latitudeR = ee.Image.pixelLonLat().select('latitude').updateMask(maskREMA);
  var latMinMaxR = latitudeR.reduceRegion({reducer: ee.Reducer.minMax(), geometry: image.geometry(), scale: 500, crs: 'EPSG:3031'});
  var longitudeR = ee.Image.pixelLonLat().select('longitude').updateMask(maskREMA);
  var lonMinMaxR = longitudeR.reduceRegion({reducer: ee.Reducer.minMax(), geometry: image.geometry(), scale: 500, crs: 'EPSG:3031'});
  var latDistR = ee.Number(latMinMaxR.get('latitude_max')).subtract(ee.Number(latMinMaxR.get('latitude_min')));
  var lonDistR = ee.Number(lonMinMaxR.get('longitude_max')).subtract(ee.Number(lonMinMaxR.get('longitude_min')));
  
  //// CS size
  var latitudeC = ee.Image.pixelLonLat().select('latitude').updateMask(maskREMA).updateMask(maskCS);
  var latMinMaxC = latitudeC.reduceRegion({reducer: ee.Reducer.minMax(), geometry: image.geometry(), scale: 500, crs: 'EPSG:3031'});
  var longitudeC = ee.Image.pixelLonLat().select('longitude').updateMask(maskREMA).updateMask(maskCS);
  var lonMinMaxC = longitudeC.reduceRegion({reducer: ee.Reducer.minMax(), geometry: image.geometry(), scale: 500, crs: 'EPSG:3031'});
  var latDistC = ee.Number(latMinMaxC.get('latitude_max')).subtract(ee.Number(latMinMaxC.get('latitude_min')));
  var lonDistC = ee.Number(lonMinMaxC.get('longitude_max')).subtract(ee.Number(lonMinMaxC.get('longitude_min')));
  
  var latFrac = latDistC.divide(latDistR);
  var lonFrac = lonDistC.divide(lonDistR);
  return image.set({latFrac:latFrac, lonFrac:lonFrac});
})};

exports.filterByCryoSatCoverage = filterByCryoSatCoverage;




var filterByREMAMosaic = function(elevMosaic,REMAMosaic){
  
  /* Mask out points where the abs. elev. diff. between REMA and REMA Mosaic is greater than *elevMosaic*
  */
  
  return(
  function(image){
  var REMAMos = REMAMosaic.updateMask(image.mask());
  var diff = REMAMos.subtract(image).abs().rename('diff');
  image = image.updateMask(diff.select('diff').lte(elevMosaic));
  return image;
})};

exports.filterByREMAMosaic = filterByREMAMosaic;



/////////////////////////////////////////////////////////////////////////////
///////////////////////////// CO-REGISTRATION ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////


var tiltBiasCorrection = function(CS2,allowedTimeSpan){
  
  /* The big funky tilt and bias correction which returns a new band to the image with the correction to be applied
  */
  
  return(
  function(image){
  var firstREMA = image;
  var imageDate = ee.Date(firstREMA.get('system:time_start'));
  var startDate = imageDate.advance(-allowedTimeSpan, 'month');
  var endDate = imageDate.advance(allowedTimeSpan, 'month');
  var CS2_image = CS2.filterDate(startDate,endDate).filterBounds(image.geometry());
  CS2_image = CS2_image.reduceToImage(['h'], ee.Reducer.mean()).rename('elevation');
  var maskREMA = firstREMA.mask();
  var maskCS = CS2_image.mask();
  firstREMA = firstREMA.updateMask(maskCS);
  // var diff = CS2_image.subtract(firstREMA).abs();
  // diff = diff.updateMask(maskREMA).updateMask(maskCS).rename('diff');
  // CS2_image = CS2_image.updateMask(diff.select('diff').lte(outlierCriteria));
  CS2_image = CS2_image.updateMask(maskREMA);
  var reg = firstREMA.geometry();//.intersection(IS);
  var dz = CS2_image.subtract(firstREMA).rename('bias');
  var longitude = ee.Image.pixelLonLat().select('longitude');
  longitude = longitude.updateMask(dz.mask());
  var latitude = ee.Image.pixelLonLat().select('latitude');
  latitude = latitude.updateMask(dz.mask());
  var constant = ee.Image.constant(1).toFloat().select('constant');
  constant = constant.updateMask(dz.mask());



  // Convert the above images into arrays, in order to construct the matrix
  var longitudeList = longitude.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100
  }).values();
  
  var latitudeList = latitude.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100
  }).values();
  
  var constantList = constant.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100
  }).values();
  

  var dzList = dz.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100
  }).values();
  
  dzList = ee.Array.cat([dzList],1).transpose(); // Making sure that the shape is nx1
  // Create matrix A
  var A = ee.Array.cat([longitudeList,latitudeList,constantList],0).transpose(); // Making sure that the shape is nx3

  // Do the matrix operation and get the coefficient vector
  var coefVec = A.matrixSolve(dzList);
  
  coefVec = coefVec.project([0]); // Go from 3x1 matrix to 3 vector
  var a = coefVec.get([0]); // x correction
  var b = coefVec.get([1]); // y correction
  var c = coefVec.get([2]); // vertical correction
  
  // Given the equation dz = ax + by + c
  
  // Calculate ax as an image
  var lon = ee.Image.pixelLonLat().select('longitude');
  lon = lon.updateMask(image.mask()).multiply(a);
  
  // Calculate bx as an image
  var lat = ee.Image.pixelLonLat().select('latitude');
  lat = lat.updateMask(image.mask()).multiply(b);
  
  // Make an image with c
  var cst = ee.Image.constant(c).toFloat().select('constant');
  cst = cst.updateMask(image.mask());
  
  // Apply the equation: dz = ax + by + c
  var z_diff = lon.add(lat).add(cst).rename('correction');
  return image.addBands(z_diff);
})};

exports.tiltBiasCorrection = tiltBiasCorrection;



var finalCorrection = function(image) {
  
  /* Apply the correction to the REMA elevations
  */
  
  var elevation = image.select('elevation');
  var corr = image.select('correction');
  return elevation.add(corr).rename('elevation')
    .set("system:time_start", image.get("system:time_start"), "latFrac", image.get("latFrac"), "lonFrac", image.get("lonFrac"), "matchups", image.get("matchups"), "percentageCS2coverage", image.get("percentageCS2coverage"));
};

exports.finalCorrection = finalCorrection;

var finalCorrectionOverlap = function(image) {
  
  /* Apply the correction to the REMA elevations
  */
  
  var elevation = image.select('elevation');
  var corr = image.select('correction');
  return elevation.add(corr).rename('elevation')
    .set("system:time_start", image.get("system:time_start"),'system:index',image.get('system:index'),'OverlapPercent',image.get('OverlapPercent'), "latFrac", image.get("latFrac"), "lonFrac", image.get("lonFrac"), "matchups", image.get("matchups"), "percentageCS2coverage", image.get("percentageCS2coverage"));
};

exports.finalCorrectionOverlap = finalCorrectionOverlap;


var tiltBiasCorrectionOverlap = function(overlapREMA){
  
  /* The big funky tilt and bias correction which returns a new band to the image with the correction to be applied
  */
  
  return(
  function(image){
  // var image = REMA.first();
  var firstREMA = image;
  
  var CS2_image = overlapREMA;
  
  var maskREMA = firstREMA.mask();
  var maskCS = CS2_image.mask();
  firstREMA = firstREMA.updateMask(maskCS);
  var diff = CS2_image.subtract(firstREMA).abs();
  diff = diff.updateMask(maskREMA).updateMask(maskCS).rename('diff');
  var reg = firstREMA.geometry();//.intersection(IS);
  var dz = CS2_image.subtract(firstREMA).rename('bias');
  var longitude = ee.Image.pixelLonLat().select('longitude');
  longitude = longitude.updateMask(dz.mask());
  var latitude = ee.Image.pixelLonLat().select('latitude');
  latitude = latitude.updateMask(dz.mask());
  var constant = ee.Image.constant(1).toFloat().select('constant');
  constant = constant.updateMask(dz.mask());
  
  // Map.addLayer(dz,{min:-4,max:4,palette:'red,white,blue'},'bias')
  
  // Convert the above images into arrays, in order to construct the matrix
  var longitudeList = longitude.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100,
    bestEffort: true
  }).values();
  
  var latitudeList = latitude.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100,
    bestEffort: true
  }).values();
  
  var constantList = constant.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100,
    bestEffort: true
  }).values();
  
  
  var dzList = dz.reduceRegion({
    reducer: ee.Reducer.toList(),
    geometry: dz.geometry(),
    crs: 'EPSG:3031',
    scale:100,
    bestEffort: true
  }).values();
  
  dzList = ee.Array.cat([dzList],1).transpose(); // Making sure that the shape is nx1
  // Create matrix A
  var A = ee.Array.cat([longitudeList,latitudeList,constantList],0).transpose(); // Making sure that the shape is nx3
  
  // Do the matrix operation and get the coefficient vector
  var coefVec = A.matrixSolve(dzList);
  
  coefVec = coefVec.project([0]); // Go from 3x1 matrix to 3 vector
  var a = coefVec.get([0]); // x correction
  var b = coefVec.get([1]); // y correction
  var c = coefVec.get([2]); // vertical correction
  
  // Given the equation dz = ax + by + c
  
  // Calculate ax as an image
  var lon = ee.Image.pixelLonLat().select('longitude');
  lon = lon.updateMask(image.mask()).multiply(a);
  
  // Calculate bx as an image
  var lat = ee.Image.pixelLonLat().select('latitude');
  lat = lat.updateMask(image.mask()).multiply(b);
  
  // Make an image with c
  var cst = ee.Image.constant(c).toFloat().select('constant');
  cst = cst.updateMask(image.mask());
  
  // Apply the equation: dz = ax + by + c
  var z_diff = lon.add(lat).add(cst).rename('correction');
  return image.addBands(z_diff);
})};

exports.tiltBiasCorrectionOverlap = tiltBiasCorrectionOverlap;


var percentOverlap = function(overlapREMA){
  return(
    function(image){
        var firstREMA = image;
  
        var CS2_image = overlapREMA;
        
        var maskREMA = firstREMA.mask();
        var maskCS = CS2_image.mask();
        
        
        // var OVdiff = CS2_image.subtract(image).abs();
        // OVdiff = OVdiff.updateMask(maskREMA).updateMask(maskCS).rename('diff');
        var reg = image.geometry();
        var scale = 500;
        // var constantImage = ee.Image.constant(1).updateMask(maskREMA).updateMask(maskCS);
        
        // Map.addLayer(constantImage)
        var areaOverlap = ee.Image.constant(1).updateMask(maskREMA).updateMask(maskCS).reduceRegion({
          reducer:ee.Reducer.sum(), 
          geometry:reg,
          scale:scale,
          crs:'EPSG:3031',
          tileScale:16})
          .get('constant');
          
        // print('area overlap', areaOverlap)
        
        var areaREMA = ee.Image.constant(1).updateMask(maskREMA).reduceRegion({
          reducer:ee.Reducer.sum(), 
          geometry:reg,
          scale:scale,
          crs:'EPSG:3031',
          tileScale:16})
          .get('constant');
          
        // print('area REMA',areaREMA)
        
        var Perc = ee.Number(areaOverlap).divide(ee.Number(areaREMA)).multiply(ee.Number(100))
    
        return image.set('OverlapPercent',Perc)
    })};
    
exports.percentOverlap = percentOverlap;
    
    
    

