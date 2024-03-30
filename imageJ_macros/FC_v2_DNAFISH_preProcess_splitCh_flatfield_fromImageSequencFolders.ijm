//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Finn Clark, Lionnet Lab, Updated 08/29/23

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~OVERVIEW~~~~~~~~~~~~~~~~

// - Damn, you accidentally saved you hyperstacks as image sequences. Time to fix it
// - This macro searches a directory for image sequences up two folders deep with n channels and then processes them. 
// - All img channels are split, flatfield corrected, and then saved to a 
// "split channels" folder within the output directory that you select. 
// - A max z projection of the processed DAPI channel is saved to a sub directory called "max_z"
// - An empty "cell mask" folder is created for downstream use

// YOU MUST SPECIFY THE CHANNEL NUMBER OF YOUR DAPI CHANNEL HERE (or whatever channel you are using to make masks).
// You must open one of your sequnces in FIJI and figure out which stack to hyperstack settings you need. Plug them in here.
// There is a chance you have to change slicing order. If so, scroll down to the run stack to hyperstack command and change it manually. DONT 
// OVERWRITE THIS MACRO
maskChan = 1;
nCh = 4;
nZSlices = 17;
nFrames = 1;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


input = getDirectory("choose INPUT directory");

//where the analysis images will go 

output = input +  "/split_channels/";

File.makeDirectory(output);

maskDir = input + "/cell_masks/";
File.makeDirectory(maskDir);

zMaxDir = input + "/max_z/"; 
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
	        			list[ctr] = input+tmplist[k]+tmplist2[l];
	        		

	        			//print(list[ctr]);
	        			ctr++;
	        			break;
        			}
				}
			}
			if((indexOf(tmplist2[l], ".tif") != -1) || (indexOf(tmplist2[l], ".lsm") != -1) || (indexOf(tmplist2[l], ".czi") != -1)){
	        	list[ctr] = input+tmplist[k];
	        	print(list[ctr]);
	        	ctr++;
	        	break
        	}
		}	
	}  

    if((indexOf(tmplist[k], ".tif") != -1) || (indexOf(tmplist[k], ".lsm") != -1) || (indexOf(tmplist[k], ".czi") != -1)){
	      list[ctr] = input;
	      //print(list[ctr]);
	      ctr++;
    }
         
}

list = Array.trim(list,ctr);

print("Number of Images Loaded:");
print(list.length);


 
//print file names in log window
print("Sequence directories to process");
for (k=0; k<list.length; k++) {
	print(list[k]);
}

//now that file names for the conditions are all combined, lets do the processing

for (i = 0; i < list.length; i++){	
    
    //open image, split channels and save each
    File.openSequence(list[i]); 
	title = getTitle();
	
	
	print("Processing");
	print(title);
	
	
	run("Stack to Hyperstack...", "order=xyzct channels="+nCh+ " slices=" + nZSlices+  " frames=" + nFrames +" display=Color");
	getDimensions(dummy,dummy,channelcount,dummy,dummy);
	print(channelcount);
	
	run("Split Channels");
	
	//saveName = File.getParent(list[i]) + i;
	
	for (myCh = 1; myCh < channelcount+1; myCh++) {
		
		if (myCh == maskChan) {
			
			selectImage("C" + maskChan + "-"  + title);
			

			// save c1
			saveName = "C1_" + i;
			saveAs("Tiff", output+ saveName);
			//save c1 for mask making
			run("Z Project...", "projection=[Max Intensity]");
			saveAs("Tiff", zMaxDir+ saveName+"MAX_Z");
			close();
		}
		else {
			
			selectImage("C" + myCh + "-" + title);
			run("Pseudo flat field correction", "blurring=50 hide stack");

			saveName = "C" + myCh + "_" + i;
			saveAs("Tiff", output+ saveName);
			close();
		}

	}
	
}




setBatchMode(false); 
        
