import qupath.lib.roi.*
import qupath.lib.objects.*
import static qupath.lib.classifiers.PathClassifierTools.*
import qupath.lib.scripting.QPEx


def starttime = new Date()
def image_name  = getCurrentImageData().getServer().getShortServerName()
// do we want to generate output? yes/no
def out = true

print "Working on image "+ image_name

def ana_list = ["PDAC-PGP01", "PDAC-PGP02", "PDAC-PGP03", "PDAC-PGP05", "PDAC-PGP06", "PDAC-PGP09", "PDAC-PGP10", "PDAC-PGP11", "PDAC-PGP13", "PDAC-PGP16", "PDAC-PGP18", "PDAC-PGP19", "PDAC-PGP20", "PDAC-PGP21", "PDAC-PGP22", "PDAC-PGP23", "PDAC-PGP24", "PDAC-PGP25", "PDAC-PGP26", "PDAC-PGP27", "PDAC-PGP28", "PDAC-PGP29", "PDAC-PGP30", "PDAC-PGP31", "PDAC-PGP32", "PDAC-PGP33", "PDAC-PGP34", "PDAC-PGP35", "PDAC-PGP36", "PDAC-PGP37", "PDAC-PGP39", "PDAC-PGP40", "PDAC-PGP47", "PDAC-PGP48", "PGP49", "PGP50", "PGP51", "PGP52", "PGP53", "PGP55", "PGP56", "PGP57", "PGP58", "PGP59", "PGP061", "PGP062", "PGP63", "PGP64", "PGP66", "PGP67", "PGP68", "PGP69", "PGP70", "PGP71", "PGP72", "PGP74", "PGP76", "PGP77", "PGP78", "PGP79", "PGP80", "PGP81", "PGP82", "PGP83", "PGP84", "PGP85", "PGP86", "PGP87", "PGP89", "PGP90", "PGP91", "PGP92", "PGP93", "PGP94", "PGP95", "PGP96", "PGP97", "PGP98", "PGP99", "PGP101", "PGP102", "PGP104", "PGP105", "PGP106", "PGP107", "PGP108", "PGP109", "PGP110", "PGP111", "PGP112", "PGP113", "PGP114", "PGP115", "PGP116", "PGP117", "PGP118", "PGP119", "PGP120", "PGP122", "PGP123", "PGP125", "PGP126", "PGP127", "PGP128", "PGP129", "PGP131", "PGP132", "PGP133", "PGP134", "PGP136", "PGP139", "PGP140", "PGP141", "PGP142", "PGP143", "PGP144", "PGP145", "PGP146", "PGP147", "PGP148", "PGP150", "PGP151", "PGP152", "PGP153", "PGP154", "PGP155", "PGP156", "PGP157", "PGP158", "PGP159", "PGP160", "PGP161", "PGP162", "PGP163", "PGP164", "PGP165", "PGP166", "PGP167", "PGP168", "PGP169", "PGP170", "PGP172", "PGP173", "PGP174", "PGP175", "PGP176", "PGP177", "PGP178", "PGP179", "PGP180", "PGP181", "PGP182", "PGP183", "PGP184", "PGP187", "PGP189", "PGP190", "PGP191", "PGP192", "PGP193", "PGP194", "PGP195", "PGP196", "PGP198", "PGP200", "PGP201", "PGP202", "PGP203", "PGP204", "PGP205", "PGP207", "PGP208", "PGP210", "PGP211", "PGP213", "PGP214", "PGP215", "PGP216", "PGP217", "PGP218", "PGP219", "PGP220", "PGP221", "PGP222", "PGP223", "PGP224", "PGP225", "PGP226", "PGP228", "PGP229"]

if(!(image_name in ana_list)){
    print "Image " + image_name + " is not in analysis list "
    print "Skipping analysis for this image"
    return
}


// List of images that we do not want to be classified
def list_skip = []

if(image_name in list_skip){
    print "Image " + image_name + " is in immune cell free list"
    print "Skipping analysis for this image"
    return
}

// check if helper regions ROI annotated on slide?
def hreg=false

def helper_annotations = getAnnotationObjects().findAll{it.getPathClass()==getPathClass("ROI")}

if(!helper_annotations.isEmpty()){
    print "Found helper anntation in "+image_name
    print "Analysis limited to helper region"
    hreg=true
}else{
    print "Did not find helper annotations named \"ROI\"."
    print "Working whole slide image."
    print "This might take a lot longer ..."
    
}


// set classifier string here
def classifier = buildFilePath(buildFilePath(PROJECT_BASE_DIR,"classifiers"),"immunecell_v2.qpclassifier")

def ignoreTheseinSubtraction=[getPathClass("totalTissue"),getPathClass("ROI"),getPathClass("Results")]

def allDetections = getDetectionObjects()

// container of detected cells
def cells = []

// run cell detection if we dont have any cells yet
if(allDetections.isEmpty()){

    //remove all objects that should only be created as helpers in this script
    def helpers = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("stromaAna") }
    
    removeObjects(helpers, true)
    print "Removed " +helpers.size()+" annotations of type StromaAna."
    
    // total Tissues dependend on whether we have defined a helper annotation or not
    def totalTissues = []
    
    if(hreg){
	
	totalTissues = helper_annotations
	
    }else{
	
	totalTissues = getAnnotationObjects().findAll {it.isAnnotation() && it.getPathClass() == getPathClass("totalTissue") && it.getROI().respondsTo('getArea')}
	
    }
    
    if(totalTissues.isEmpty()){
	print "Did not find an annotation to analyse"
	return
    }

    subtractions = getAnnotationObjects().findAll{it.isAnnotation() && !(it.getPathClass() in ignoreTheseinSubtraction) && it.getROI().respondsTo('getArea')}
    
    for (totalTissue in totalTissues){
	def total = []
	def polygons = []
	
	for (subtractyBit in subtractions){
	    
	    // This is conversion from Area to polygon
	    if (subtractyBit instanceof AreaROI){
		subtractionROIs = PathROIToolsAwt.splitAreaToPolygons(subtractyBit.getROI())
		total.addAll(subtractionROIs[1])
	    } else {total.addAll(subtractyBit.getROI())}
	}
	// Conversion from Area to polygon
	if (totalTissue instanceof AreaROI){
	    polygons = PathROIToolsAwt.splitAreaToPolygons(totalTissue.getROI())
	    total.addAll(polygons[0])
	} else { polygons[1] = totalTissue.getROI()}
	
	def newPolygons = polygons[1].collect {
	    updated = it
	    for (hole in total)
	    updated = PathROIToolsAwt.combineROIs(updated, hole, PathROIToolsAwt.CombineOp.SUBTRACT)
	    return updated
	}
	
	// add new ones
	annotations = newPolygons.collect {new PathAnnotationObject(updated,getPathClass("stromaAna"))}
	
	addObjects(annotations)
    }
 
    
    selectObjects{ it.isAnnotation() && it.getPathClass() == getPathClass("stromaAna")};
    
    runPlugin('qupath.imagej.detect.nuclei.WatershedCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 0.1,  "maxBackground": 2.0,  "watershedPostProcess": true,  "excludeDAB": false,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

    //Store the annotations in "annotations" and create an empty array where you can keep the list of cells you want to classify.
    hierarchy = getCurrentHierarchy();
    annotationsA = hierarchy.getSelectionModel().getSelectedObjects()
    //Iterate over all annotations and keep track of the cells within them in the "cells" list
    
    for (annotation in annotationsA){
	cells.addAll(hierarchy.getDescendantObjects(annotation, null, PathDetectionObject));
    }
}
// alternatively analyse all existing detections
else{
    print "Found detections on slide."
    print "Analyse existing detections ..."
    cells=getDetectionObjects();
}


// calculate the area we are examining
// Area of analysed ROI
def totalArea = 0.0

stromaAnas = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("stromaAna")}


for(area in stromaAnas){
    
    totalArea = totalArea + area.getROI().getArea()
}


// Apply classifier to all detections inside the selected annotations
// run cell detection on stromaAna
c1 = loadClassifier( new File(classifier) )

c1.classifyPathObjects(cells);


// only add cells if they are not added yet
if(allDetections.isEmpty()){
    addObjects(cells)
}

//write statistics to file
if(out){    
// prepare output
    def outdir_tumor="tumorData"
    def outdir_immune="immuneData"
        
    def foutdir_tumor = buildFilePath(PROJECT_BASE_DIR,outdir_tumour)
    def foutdir_immune = buildFilePath(PROJECT_BASE_DIR,outdir_immune)
    if(!QPEx.fileExists(foutdir_tumor)){
	QPEx.mkdirs(foutdir_tumor)
    }
    if(!QPEx.fileExists(foutdir_immune)){
	QPEx.mkdirs(foutdir_immune)
    }
        
    def str_filename  = getCurrentImageData().getServer().getShortServerName()+".data"
    def filename_tumor = buildFilePath(foutdir_tumor,str_filename)
    def filename_immune = buildFilePath(foutdir_immune,str_filename)
    
    fileT=new File(filename_tumor)
    wT = fileT.newWriter()
    fileI=new File(filename_immune)
    wI = fileI.newWriter()
}

def immuneCells = []

for(cell in cells){
    if(cell.getPathClass() == getPathClass("ImmuneCells")){
	immuneCells.add(cell)
    }
}

print "Starting to analyse the tumor immune cell distance"

annos = getAnnotationObjects().findAll{it.isAnnotation() && it.getROI().respondsTo('getArea')}

for immuneCell in immuneCells{
	if(immuneCell.respondsTo("getNucleusROI") && immuneCell.getNucleusROI() != null){
                def roi = immuneCell.getNucleusROI()
                double cx = roi.getCentroidX().round(1)
	        double cy = roi.getCentroidY().round(1)
	        if(roi.getCentroidX() instanceof Number && roi.getCentroidY() instanceof Number){
		        wI << cx<< ", "<< cy <<"\n"
	        }
        }
}

for (anno in annos){
    def roi = anno.getROI()
    
    if(roi == null){
	return
    }
    
    points = roi.getPolygonPoints()
    for (point in points) {
        p_x = point.getX()
        p_y = point.getY()
        point_string = p_x + ", " + p_y
        wT << point_string << System.lineSeparator()
    }                     
    wT << '#' << anno.getPathClass() << System.lineSeparator()
}


// if output was generated
if(out){
    wT.close()
    wI.close()
}
def endtime = new Date()
def minutes = (endtime.getTime()-starttime.getTime())/60.0/1000.0

print ''
print('Done image '+image_name+' in ' + minutes + 'minutes.')
print ''
