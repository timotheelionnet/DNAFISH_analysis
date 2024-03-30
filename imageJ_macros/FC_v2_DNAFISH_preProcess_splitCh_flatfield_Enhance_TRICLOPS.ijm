//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Finn Clark, Lionnet Lab, Updated 12/30/22

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~OVERVIEW~~~~~~~~~~~~~~~~

// - This macro searches a directory for hyperstacks up two folders deep with n channels and then processes them. 
// - All img channels are split, flatfield corrected, and then saved to a 
// "split channels" folder within the output directory that you select. 
// - A max z projection of the processed DAPI channel is saved to a sub directory called "max_z"
// - An empty "cell mask" folder is created for downstream use

// YOU MUST SPECIFY THE CHANNEL NUMBER OF YOUR DAPI CHANNEL HERE (or whatever channel you are using to make masks).
maskChan = 1;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


input = getDirectory("choose INPUT directory");

//where the analysis images will go 
outputRoot = File.getParent(input);
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

print("Number of Images Loaded:");
print(list.length);


 
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
		
		print("Processing");
		print(title);
		
		run("Split Channels");
		
		for (myCh = 1; myCh < channelcount+1; myCh++) {
			
			if (myCh == maskChan) {
				
				selectImage("C" + maskChan + "-"  + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				

				// save c1
				tmpimg = getTitle();
				saveAs("Tiff", output+tmpimg);
				//save c1 for mask making
				run("Z Project...", "projection=[Max Intensity]");
				saveAs("Tiff", zMaxDir+tmpimg+"MAX_Z");
				close();
			}
			else {
				
				selectImage("C" + myCh + "-" + title);
				run("Pseudo flat field correction", "blurring=50 hide stack");
				run("Enhance Contrast...", "saturated=0.001 process_all use");	

				tmpimg = getTitle();
				saveAs("Tiff", output+tmpimg);
				close();
			}

		}
	}
}




setBatchMode(false); 
        
