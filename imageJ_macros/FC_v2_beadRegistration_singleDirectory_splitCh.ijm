//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Finn Clark, Lionnet Lab, Updated 12/30/22

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~OVERVIEW~~~~~~~~~~~~~~~~

// - This macro is for processing a bead scan acquisition before downstream analysis.
// - The input must be a folder containing the output stacks from a bead scan.
// - The macro applies a flat field correction, then normalizes, then saves split channels in a separate folder.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nChannels = 3;
input = getDirectory("choose your directory");

//where the analysis images will go 
output = input +"/split_channels/";
File.makeDirectory(output);
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
 
// I want to go through and do stack to hyperstack for all images
//nChannels = 3; 
//newNSlices = nSlices/nChannels;
//
//run("Stack to Hyperstack...", "order=xyczt(default) channels=" +nChannels + "slices=" newNSlices+ "frames=1 display=Color");
//




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
		getDimensions(dummy,dummy,channelcount,dummy,dummy);
		
		run("Split Channels");
		
		for(myCh = 1; myCh < channelcount +1 ;  myCh++) {
			
			selectImage("C" + myCh +"-" + title);
			run("Pseudo flat field correction", "blurring=50 hide stack");
			tmpimg = "C" + myCh+ "_" + fov;
			saveAs("Tiff", output+tmpimg);
			close();
			
		}
    }
}

print("Operation Complete");

setBatchMode(false); 
        