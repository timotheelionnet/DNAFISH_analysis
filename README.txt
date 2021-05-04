A series of scripts for DNA FISH analysis, co-developed by Yi Fu and Timothee Lionnet  
  
The master script is DNA_FISH_pipeline.m, which exacute the following:  
generate_test_img.m  
Generate test images using functions in qc.  
  
DNA_FISH_pipeline_1.m  
Input raw DNA FISH images, perform Airlocalize (including clean up double detections) and crop loci from images.  
  
match_DNA_FISH_spots_between_channels.m  
Input cell mask images and .loc3 of centers of cropped loci, match up loci between different FISH channels according to cell ID and distances.  
  
analysis_overlap_between_loci.m  
Input loci match results and cropped loci images, compute cross correlation between match loci pairs.  
  
mensure_cropped_loci.m  
Input cropped loci images, compute centroid locations, radio of gyrations, and volumns of loci.
  
plotPDF02032021.m  
Input data, plot PDF and CDF.  
  
analysis_cell_cycle.m  
Input loci match results, identify cell cycle(G1/G2) by number of DNA FISH spots, then plot DAPI/EU/IF intensity and distances between loci center grouped by cell cycle.  
    
*fileset, useful, and AIIRLOCALIZE1_stable is too big for uploading
