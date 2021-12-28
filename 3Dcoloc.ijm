/*

#  BIOIMAGING - INEB/i3S
Eduardo Conde-Sousa (econdesousa@gmail.com)

## 3D coloc
 
### code version
1.0

### last modification
09/12/2021

### Requirements
	* CLIJ2
	* IJPB-plugins
	* 3D ImageJ Suite


### Attribution:
If you use this macro please add in the acknowledgements of your papers and/or thesis (MSc and PhD) the reference to Bioimaging and the project PPBI-POCI-01-0145-FEDER-022122.
As a suggestion you may use the following sentence:
 * The authors acknowledge the support of the i3S Scientific Platform Bioimaging, member of the national infrastructure PPBI - Portuguese Platform of Bioimaging (PPBI-POCI-01-0145-FEDER-022122).

*/

#@ File (label = "Input directory", style = "directory") inDir
#@ String (label = "File suffix", value = ".ics") EXT 
#@ Integer (label = "ring size [pixels]",description="values that correspond to more than 2 microns may result in removinginner regions", value = 2 ) borderSize 
#@ Boolean (label =  "set batch mode?", value = true) batchMode
// number of pixels on edges (1/2 in z direction)
#@ Integer (label = "Signal Channel", value = 1 ) signalChannel
#@ Integer (label = "Mitochondria Channel", value = 2 ) mitoChanel

var width;
var height;
var depth;
var unit;
var filename = "";
var outDir = "";
var is2pull = !batchMode;
var mitoImage = "";
var signalImage = "";
var Label = newArray();
var Mean = newArray();
var Std = newArray();
var	Max = newArray();
var	Min = newArray();
var	Median = newArray();
var	Mode = newArray();
var	Volume = newArray();
var	Region = newArray();
var	FileName = newArray();
	


/*
Initial setup
*/

// Init GPU
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();
// setup image variables
close("*")

if (!endsWith(inDir, File.separator )) {
	inDir = inDir + File.separator;
}
outDir = inDir + "Results" + File.separator ;
if (!File.exists(outDir)) {
	File.makeDirectory(outDir);
}
print("\\Clear");
selectWindow("Log");
filelist = getFileList(inDir) ;
for (i = 0; i < lengthOf(filelist); i++) {
    if (endsWith(filelist[i], EXT)) { 
        filename = filelist[i];
        print(filename);
        proc1file();
    } 
}

Array.show("Results", FileName,Label,Region,Mean,Std,Max,Min,Median,Mode,Volume);
saveAs("Results", outDir + File.getNameWithoutExtension(inDir) + "_consolidatedResults.tsv" );
setBatchMode(false);

function proc1file(){
	close("*");
	selectWindow("Log");
	showStatus("opening " + filename);
	showProgress(0);
	setBatchMode(batchMode);	
	run("Bio-Formats Importer", "open=["+inDir + filename+"] autoscale color_mode=Colorized rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	getVoxelSize(width, height, depth, unit);
	
	mainID = getImageID();	
	selectImage(mainID);
	run("Duplicate...", "title=mito duplicate channels="+mitoChanel);
	mitoImage = getTitle();
	selectImage(mainID);
	run("Duplicate...", "title=signal duplicate channels="+signalChannel);
	signalImage = getTitle();
	showProgress(0.1);
	
	Ext.CLIJ2_pushCurrentZStack(mitoImage);
	
	showStatus("segment mitochodria");
	// Top Hat Box
	radiusX = 3.0;
	radiusY = 3.0;
	radiusZ = 1.5;
	Ext.CLIJ2_topHatBox(mitoImage, image_1, radiusX, radiusY, radiusZ);
	Ext.CLIJ2_release(mitoImage);
	//Ext.CLIJ2_pull(image_1);
	showProgress(0.2);
	
	// Automatic Threshold
	showStatus("top hat filter");
	method = "Triangle";
	mitoMask = mitoImage+"_mask";
	Ext.CLIJ2_automaticThreshold(image_1, mitoMask, method);
	Ext.CLIJ2_saveAsTIF(mitoMask, outDir + replace(filename,EXT,"_mitoMask.tif"));
	if (is2pull){
		Ext.CLIJ2_pullBinary(mitoMask);
		setVoxelSize(width, height, depth, unit);
	}
	Ext.CLIJ2_release(image_1);
	showProgress(0.3);

	showStatus("segment signal channel");
	Ext.CLIJ2_pushCurrentZStack(signalImage);
	// Automatic Threshold
	method = "Triangle";
	signalMask= signalImage+"_mask";
	Ext.CLIJ2_automaticThreshold(signalImage,tmp, method);
	Ext.CLIJ2_release(signalImage);
	Ext.CLIJ2_maximum3DBox(tmp, signalMask, 3, 3, 1.5);
	Ext.CLIJ2_release(tmp);
	Ext.CLIJ2_saveAsTIF(signalMask, outDir + replace(filename,EXT,"_signalMaskRestricted.tif"));
	if (is2pull){
		Ext.CLIJ2_pullBinary(signalMask);
		setVoxelSize(width, height, depth, unit);
	}
	showProgress(0.4);

	
	// Binary And
	showStatus("get mask");
	maskRestricted = "maskRestricted";
	Ext.CLIJ2_binaryAnd(mitoMask, signalMask, maskRestricted);
	Ext.CLIJ2_saveAsTIF(maskRestricted, outDir + replace(filename,EXT,"_mitoMaskRestricted.tif"));
	if (is2pull){
		Ext.CLIJ2_pullBinary(maskRestricted);
		setVoxelSize(width, height, depth, unit);
	}
	showProgress(0.5);
	
	// get final masks
	mitoCentralMask = "mitoCentralMask";
	mitoBoundaryMask = "mitoBoundaryMask";
	mitoOuterBoundaryMask = "mitoOuterBoundaryMask";
	
	showProgress(0.6);
	Ext.CLIJ2_minimum3DSphere(maskRestricted, mitoCentralMask, borderSize, borderSize, borderSize/2);
	Ext.CLIJ2_maximum3DSphere(maskRestricted, image_2, borderSize, borderSize, borderSize/2);
	
	
	Ext.CLIJ2_binarySubtract(maskRestricted, mitoCentralMask, mitoBoundaryMask);
	Ext.CLIJ2_binarySubtract(image_2,maskRestricted, mitoOuterBoundaryMask);
	Ext.CLIJ2_release(image_2);
	
	if (is2pull){
		Ext.CLIJ2_pullBinary(mitoBoundaryMask);
		setVoxelSize(width, height, depth, unit);
		Ext.CLIJ2_pullBinary(mitoCentralMask);
		setVoxelSize(width, height, depth, unit);
		Ext.CLIJ2_pullBinary(mitoOuterBoundaryMask);
		setVoxelSize(width, height, depth, unit);
	}
	showProgress(0.7);
	
	mitoRegions = "mitoRegions";
	Ext.CLIJ2_addImagesWeighted(mitoCentralMask, mitoBoundaryMask, tmp, 1, 2);
	Ext.CLIJ2_addImagesWeighted(tmp, mitoOuterBoundaryMask, mitoRegions, 1, 3);
	Ext.CLIJ2_pull(mitoRegions);
	setMinAndMax(0, 3);
	run("glasbey_on_dark");
	
	Ext.CLIJ2_release(mitoMask);
	Ext.CLIJ2_release(signalMask);
	Ext.CLIJ2_release(maskRestricted);
	Ext.CLIJ2_clear();
	
	showProgress(0.8);
	run("Tile");
	
	showStatus("getting results");
	run("Intensity Measurements 2D/3D", "input=signal labels=mitoRegions mean stddev max min median mode volume");
	Table.setColumn("Region", newArray("Center","Boundary","Outside"));
	showProgress(0.9);
	t 		= Table.getColumn("Label");
	Label 	= Array.concat(Label ,t);
	t 		= Table.getColumn("Mean");
	Mean	= Array.concat(Mean ,t);
	t		= Table.getColumn("StdDev");
	Std		= Array.concat(Std ,t);
	t		= Table.getColumn("Max");
	Max		= Array.concat(Max ,t);
	t		= Table.getColumn("Min");
	Min		= Array.concat(Min ,t);
	t		= Table.getColumn("Median");
	Median	= Array.concat(Median ,t);
	t		= Table.getColumn("Mode");
	Mode 	= Array.concat(Mode ,t);
	t 		= Table.getColumn("Volume");
	Volume 	= Array.concat(Volume ,t);
	t		= Table.getColumn("Region");
	Region	= Array.concat(Region ,t);
	FileName=Array.concat(FileName ,newArray(filename,filename,filename));	
	saveAs("Results", outDir + replace(filename,EXT,".tsv"));
	showProgress(1);
	selectWindow(replace(filename,EXT,".tsv"));
	run("Close");
	close("*");
}
