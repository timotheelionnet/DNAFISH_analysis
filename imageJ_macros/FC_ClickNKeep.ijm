// Finn Clark 3/13/2024
// You can click on cells using multi poitn tool
// Then you can save the table for downstream filtering 
// Open all images of interest
// Use multi point to mark cells of interest
// Run this macro once per FOV
// Sometimes it will ask if you want to save changes, no need to

// we will store our results in the parent dir of the file
cur_dir = getDirectory("image");
fov_name = getTitle();

fov_name_no_ext = fov_name.substring(0, fov_name.lastIndexOf("."));


//meaures points and save as csv
run("Measure");
save_name = fov_name_no_ext + "_clickNkeep.csv";
saveAs("Results", cur_dir + save_name);

// close the current image
close();

// close the results window
selectWindow("Results");

run("Close");