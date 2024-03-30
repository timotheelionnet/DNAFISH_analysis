//Batch Processing: Z projection of max intensity in stack or hyperstack

//root folder where all data is
input = getDirectory("choose your directory");
//input = "/Users/lionnett/Dropbox (HHMI)/Training Set Data (0831-0903)/";

//where the analysis images will go 
output = input + "/output/";
File.makeDirectory(output);
//output = "/Users/lionnett/Documents/junk/";
//batch mode so nothing is displayed to accelerate execution
setBatchMode(true); 

//create analysis directories if they do not exist

//collect file and folder names from root dir
tmplist = getFileList(input);

//declare arbitrarily large file list
list = newArray(10000);

//counter of the image files found in the directory and subdirectories
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
 
 print("list of files to treat:");
 for (k=0; k<list.length; k++) {
 	print(list[k]);
 }
 
print("starting file processing...");
//now that file names for the conditions are all combined, lets do the analysis
for (i = 0; i < list.length; i++){
		//print file name in log window
		//print(d2s(indexOf(list[i], ".tif"),1));
		
        if ((indexOf(list[i], ".tif") != -1) || (indexOf(list[i], ".lsm") != -1) || (indexOf(list[i], ".czi") != -1)) { // making sure it is an image
	        //open image, split channels and save each
	        open(list[i]); 
        	title = getTitle();
        	extIndex = indexOf(title, ".");
        	print("processing file "+ title);
			getDimensions(dummy,dummy,channelcount,dummy,dummy);

			run("Z Project...", " projection=[Max Intensity]");

			saveTitle = substring(title, 0, extIndex);
			saveAs("Tiff", output+saveTitle+"MAX_Z");
			close();

				
			}
			else {
				close();
			}
			close("*");
		}
close("*");

setBatchMode(false); 
        