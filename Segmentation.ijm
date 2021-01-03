close("*");
run("Clear Results");
inputDirectory = getDirectory("Choose a Directory of Images");
fluoImageDirectory = inputDirectory + "fluo-4" ;
//open(fluoImageDirectory+"/output/Average.tif");
filelist = getFileList(fluoImageDirectory);
setBatchMode(true);
open(inputDirectory + File.separator + "fluo-4" + File.separator + filelist[0]);
rename("image1")
for (i=1; i< filelist.length; i++){
	if(endsWith(filelist[i], ".tif")){
		open(inputDirectory + File.separator + "fluo-4" + File.separator + filelist[i]);
		rename("image2");
		imageCalculator("Average create", "image1","image2");;
		rename("image1");
		close("\\Others");
	}
}
setBatchMode(false);
f = File.open(inputDirectory + "Results.txt");
rename("org");
save(inputDirectory + "\\Average.tif");
run("Duplicate...", " ");
rename("WholeCell");
run("Enhance Contrast...", "saturated=0.3 normalize");
run("Gaussian Blur...", "sigma=2");
run("High pass", "radius=50 shift");
run("Subtract Background...", "rolling=50");
//run("Threshold...");
setAutoThreshold("Huang dark");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Median...", "radius=5");
save(inputDirectory+"\\WholeCellSeg.tif");


filelist = getFileList(inputDirectory + File.separator + "nuclei");
for (i=0; i< filelist.length; i++){
	if(endsWith(filelist[i], ".tif")){
		open(inputDirectory + File.separator + "nuclei" + File.separator + filelist[i]);
		rename("nuclei");
		break;
	}
}
run("Enhance Contrast...", "saturated=0.3 normalize");
run("Gaussian Blur...", "sigma=2");
run("High pass", "radius=50 shift");
run("Subtract Background...", "rolling=50");
//run("Threshold...");
setAutoThreshold("Otsu dark");
setOption("BlackBackground", true);
run("Convert to Mask");
run("Median...", "radius=5");
run("Watershed");
run("Set Measurements...", "area mean centroid perimeter redirect=None decimal=3");
run("Analyze Particles...", "size=500-Infinity display clear add");
selectWindow("WholeCell");
run("Duplicate...", " ");
rename("contour");
run("Find Edges");

run("Marker-controlled Watershed", "input=contour marker=nuclei mask=WholeCell binary calculate use");
setOption("ScaleConversions", true);
run("8-bit");
setThreshold(1, 255);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Set Measurements...", "area mean centroid perimeter area_fraction redirect=nuclei decimal=3");
run("Analyze Particles...", "size=500-Infinity display clear add");

a = nResults;		
roiManager("show all with labels");							
roiManager("multi-measure measure_all append");
//roiManager("save", inputDirectory+"\\ROI.zip");
for(i=0; i<a; i++)
	print(f, (i+1)+"\t"+ getResult("Area", i)+"\t"+ 
//		getResult("Area", i+2*a)*getResult("%Area", i+2*a)/100+"\t"+
	getResult("%Area", i)/100+"\t"+
//		getResult("%Area", i+2*a)/100+"\n");
			"\n");
print(f, "\n");
File.close(f);

roiManager("show none");
roiManager("deselect");
Nuc = newArray(a);
ignoredList = newArray(a);
Array.fill(Nuc, -1);
Array.fill(ignoredList, -1);
counterNuc = 0;
counterIgnoredList = 0;
MinAreaR = 0.1;
MinArea = 800;
for(i=0; i<a; i++){ 			
	if(getResult("%Area", i)/100 >= MinAreaR && getResult("Area", i)>=MinArea){
		print("Nuc "+getResult("%Area", i)/100+" "+getResult("Area", i));
		Nuc[i]= i;
	// 				setColor("green");
	// 				drawString(getResult("%Area", i)/100 
	// 	//						+"\n"+getResult("%Area", i+a)/100
	// 							+"\n"+getResult("Area", i), 
	// 							getResult("BX", i), getResult("BY", i)+getResult("Height", i)/2);
		counterNuc++;
	}
	else{
		print("ignored "+getResult("%Area", i)/100+" " +getResult("Area", i)/100);
		ignoredList[i] = i;
// 				setColor("red");
// 				drawString(getResult("%Area", i)/100 
//// 							+"\n"+getResult("%Area", i+a)/100
// 							+"\n"+getResult("Area", i), 
//				getResult("BX", i), getResult("BY", i)+getResult("Height", i)/2);
		counterIgnoredList++;
	}
}
Nuc = Array.deleteValue(Nuc, -1);
ignoredList = Array.deleteValue(ignoredList, -1);
roiManager("select", ignoredList);
roiManager("delete");
roiManager("save", inputDirectory + File.separator +"roi.zip");
// 		for(i=0; i<counterNuc; i++){
//	 		roiManager("select", Nuc[i]);
//	 		setForegroundColor(0, 255, 0);
//	 		run("Draw", "stack");
// 		}

selectWindow("org");
roiManager("show all with labels");

run("Enhance Contrast...", "saturated=0.3 normalize");

n = roiManager('count');
//f = File.open("/home/bill/BCH/Link to University/Harvard/iPSC-CM transients/Mike_/MT_062820_PKP2 R413X_DSC2_Y332X_ISO/scan/Well__D_003/Output.txt");
f = File.open(inputDirectory + "Output.txt");
text ='';
for(i=0; i<n; i++)
	text += ""+i+",";
text +="\n"
	print(f, text);
filelist = getFileList(fluoImageDirectory);
setBatchMode(true);
close("*");
index = 0;
i=0
while (index<filelist.length) {
	if (endsWith(filelist[index], ".tif")) { 
        open(fluoImageDirectory + File.separator + filelist[index]);
        rename(""+(i)+".tif");
        print(i+" "+fluoImageDirectory + File.separator  + filelist[index]);
		run("Set Measurements...", "area mean centroid perimeter area_fraction redirect=" +i+".tif" + " decimal=8");
		text = "";
		for (j = 0; j < n; j++) {
		    roiManager('select', j);
		    run("Measure");
		   // print(getTitle()+" "+j+" "+ getResult("Mean", nResults-1));
		    text += ""+getResult("Mean", nResults-1) + ",";
		}
		text += "\n";
		print(f, text);
		run("Clear Results");
	    close("\\Others");
	    i++;
    }   
    index++;
}
