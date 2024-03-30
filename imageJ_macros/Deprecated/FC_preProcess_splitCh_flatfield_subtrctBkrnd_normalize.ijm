//This searches a directory for hyperstacks with 3 channels, including DAPI. 
// All img channels are split, flatfield corrected, bkrnd subtracted, and normalized and then saved to subdir
//A max z projection of the processed DAPI channel is saved to a sub directory


input = getDirectory("choose INPUT directory");

//where the analysis images will go 
outputRoot = getDirectory("choose OUTPUT directory");
output = outputRoot +  "/split_channels/";

File.makeDirectory(output);

maskDir = outputRoot + "/cell_masks/";
File.makeDirectory(maskDir);

zMaxDir = outputRoot + "/max_z/"; 
File.makeDirectory(zMaxDir);
//******************************************************************************

//batch mode so nothing is displayed to accelerate execution
setBatchMode(true); 

//create analysis directories if they do not exist

tmplist = getFileList(input);

//declare arbitrary large list that will be filled with file names 
list = newArray(1000);

//counter of the file names detected
ctr = 0;


//search through file tree up to two subfolders deep for .tif images and makes file path list
for (k=0; k<tmplist.length; k++) {
 	
    if (endsWith(tmplist[k], "/") || endsWith(tmplist[k], "\\") ){
		tmplist2 = getFileList(input+tmplist[k]);
		print(tmplist[k]);	
		for (l=0; l<tmplist2.length; l++) {
			print(tmplist2[l]);	
			if (endsWith(tmplist2[l], "/") || endsWith(tmplist2[l], "\\") ){
				tmplist3 = getFileList(input+tmplist[k]+tmplist2[l]);
				for (m=0; m<tmplist3.length; m++) {	
					if( (indexOf(tmplist3[m], ".tif") != -1) || (indexOf(tmplist3[m], ".lsm") != -1) || (indexOf(tmplist3[m], ".czi") != -1)){
	        			list[ctr] = input+tmplist[k]+tmplist2[l]+tmplist3[m];
	        			//print(list[ctr]);
	        			ctr++;
        			}
				}
			}
			if((indexOf(tmplist2[l], ".tif") != -1) || (indexOf(tmplist2[l], ".lsm") != -1) || (indexOf(tmplist2[l], ".czi") != -1)){
	        	list[ctr] = input+tmplist[k]+tmplist2[l];
	        	//print(list[ctr]);
	        	ctr++;
        	}
		}	
	}  

    if((indexOf(tmplist[k], ".tif") != -1) || (indexOf(tmplist[k], ".lsm") != -1) || (indexOf(tmplist[k], ".czi") != -1)){
	      list[ctr] = input+tmplist[k];
	      //print(list[ctr]);
	      ctr++;
    }
         
 }
 list = Array.trim(list,ctr);
 
//print file names in log window
 for (k=0; k<list.length; k++) {
 	print(list[k]);
 }

//now that file names for the conditions are all combined, lets do the processing

for (i = 0; i < list.length; i++){	
        if (indexOf(list[i], ".tif") != -1) { // making sure it is an image
	        //open image, split channels and save each
	        open(list[i]); 
        	title = getTitle();
			getDimensions(dummy,dummy,channelcount,dummy,dummy);
			if (channelcount>1){
				run("Split Channels");
				selectImage("C1-" + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				run("Subtract Background...", "rolling=800 sliding stack");
				run("Enhance Contrast...", "saturated=0.00001 process_all use");
				// save c1
				tmpimg = getTitle();
				saveAs("Tiff", output+tmpimg);
				//save c1 for mask making
				run("Z Project...", "projection=[Max Intensity]");
				saveAs("Tiff", zMaxDir+tmpimg+"MAX_Z");
				close();
			
				selectImage("C2-" + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				run("Subtract Background...", "rolling=800 sliding stack");
				run("Enhance Contrast...", "saturated=0.00001 process_all use");
				tmpimg = getTitle();
				saveAs("Tiff", output+tmpimg);
				close();
				
				if (channelcount>2){

				selectImage("C3-" + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				run("Subtract Background...","rolling=800 sliding stack");
				run("Enhance Contrast...", "saturated=0.00001 process_all use");
				tmpimg = getTitle();
				saveAs("Tiff", output+tmpimg);
				close();
				}

				if (channelcount>3){

				selectImage("C4-" + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				run("Subtract Background...","rolling=800 sliding stack");
				run("Enhance Contrast...", "saturated=0.00001 process_all use");
				tmpimg = getTitle();
				saveAs("Tiff", output+tmpimg);
				close();

				}

				if (channelcount>4){

				//selectImage("C5-" + title);
				//tmpimg = getTitle();
				//saveAs("Tiff", output+tmpimg);
				//close();
				}
				
				if (channelcount>5){

				//selectImage("C6-" + title);
				//tmpimg = getTitle();
				//saveAs("Tiff", output+tmpimg);
				//close();
				}

				if (channelcount>6){

				//selectImage("C7-" + title);
				//tmpimg = getTitle();
				//saveAs("Tiff", output+tmpimg);
				//close();
				}
				
				if (channelcount>7){

				//selectImage("C8-" + title);
				//tmpimg = getTitle();
				//saveAs("Tiff", output+tmpimg);
				//close();
				}
				
				if (channelcount>8){

				//selectImage("C9-" + title);
				//tmpimg = getTitle();
				//saveAs("Tiff", output+tmpimg);
				//close();
				}

				if (channelcount>9){

				//selectImage("C10-" + title);
				//tmpimg = getTitle();
				//saveAs("Tiff", output+tmpimg);
				//close();
				}
			}
		}
}



setBatchMode(false); 
        
