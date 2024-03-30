//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Finn Clark, Lionnet Lab, Updated 12/30/22

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~OVERVIEW~~~~~~~~~~~~~~~~

// - This macro is for batch processing of bead scans. 
// - The input must be a folder containing at least one folder of stacks from a bead scan
// - The macro goes through each bead scan folder, applies a flat field correction, then normalizes, then saves split channels in that folder.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nChannels = 4;
input = getDirectory("choose your input directory"); // must be a directory containing any number of directories of bead scans 


output = input;

//******************************************************************************

//batch mode so nothing is displayed to accelerate execution
setBatchMode(true); 

//create analysis directories if they do not exist

tmplist = getFileList(input);

//declare arbitrary large list that will be filled with file names 
list = newArray(10000);
saveDirList = newArray(10000);

//counter of the file names detected
ctr = 0;


//search through file tree up to two subfolders deep for .tif images and makes file path list
for (k=0; k<tmplist.length; k++) {
 	
    if (endsWith(tmplist[k], "/") || endsWith(tmplist[k], "\\") ){
    	output = input + tmplist[k] + "split_channels/";
    	File.makeDirectory(output); 
    	
    	
		tmplist2 = getFileList(input+tmplist[k]);
 //these will be the images in the subdirectory
		print("subDir List");
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
				saveDirList[ctr] = input+tmplist[k] + "split_channels/"; //add the corresponding output subdir to list. index ctr will map image to subdir.
	        	list[ctr] = input+tmplist[k]+tmplist2[l];
//add the image file path to the list
	        	
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
 	print("here 1");
 	print(list[k]);
 	print("here 2");
 	print(saveDirList[k]);
}



//now that file names for the conditions are all combined, lets do the analysis
//for (i = 0; i < list.length; i++){
for (i = 0; i < list.length; i++){	
        if (indexOf(list[i], ".tif") != -1) { // making sure it is an image
	        //open image, run stack to hyperstack, split channels and save each
	        open(list[i]); 

			
			newNSlices = (nSlices / nChannels);

			run("Stack to Hyperstack...", "order=xyczt(default) channels=" + nChannels + " slices=" + newNSlices+ " frames=1 display=Color");
	        
        	title = getTitle();
        	fov = "FOV" + i; 
        	output = saveDirList[i];
			getDimensions(dummy,dummy,channelcount,dummy,dummy);
			
			run("Split Channels");
				
			for (myCh = 1; myCh < channelcount+1; myCh++) {
				
				
				selectImage("C" + myCh + "-" + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				
				run("Enhance Contrast...", "saturated=0.00001 process_all use");
				tmpimg = "C" + myCh + "_" + fov;;
				saveAs("Tiff", output+tmpimg);
				close();
		

			}
		
		}
}

print("Operation Complete");

setBatchMode(false); 
