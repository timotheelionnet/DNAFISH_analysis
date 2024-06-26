function analyze_multiple_images


dir = 'D:\fish\2010-03-30\Actb800-1200andExonCy3_IntronCy35_MS2Cy5\test\';


nst1 = 'st3reg.tif';
nst2 = 'st35reg.tif';
nst3 = 'st5reg.tif';

for i = 1:size(dir,1)
    %register channels
    register_mp_and_stack(dir(i,:));
    
    %% channel 1 (Cy3)
    st1 = timtiffread([dir(i,:),nst1]);
    
    p = load_parameters_for_detection2(64,64,200,2,'3D mask full',107,370,...
    2,1,1.5,1.5,3,0.01,0,'SD',20000,'D:/temp/xxx.mat',200);
    res3_3 = detection2(st1,p);%include non user interfaced option
    
    p = load_parameters_for_detection2(64,64,200,2,'3D mask full',107,370,...
    2,1,1.5,1.5,0.5,0.01,0,'SD',20000,'D:/temp/xxx.mat',200);
    res3_05 = detection2(st1,p);%include non user interfaced option
    
    save([dir(i,:),'cy3.mat'], 'res3_3', 'res3_05');
    clear('st1','res3_05','p');
    
    %% channel 2 (Cy3.5)
    st2 = timtiffread([dir(i,:),nst2]);
    
    p = load_parameters_for_detection2(64,64,200,2,'3D mask full',138,320,...
    2,1,1.5,1.5,6,0.01,0,'SD',20000,'D:/temp/xxx.mat',200);
    res35_6 = detection2(st2,p);%include non user interfaced option
    
    save([dir(i,:),'cy35.mat'], 'res35_6');
    clear('st2','p');
    
    %% channel 3 (Cy5)
    st3 = timtiffread([dir(i,:),nst3]);
    
    p = load_parameters_for_detection2(64,64,200,2,'3D mask full',128,530,...
    2,1,1.5,1.5,3,0.01,0,'SD',20000,'D:/temp/xxx.mat',200);
    res5_3 = detection2(st3,p);%include non user interfaced option
    
    p = load_parameters_for_detection2(64,64,200,2,'3D mask full',128,530,...
    2,1,1.5,1.5,0.5,0.01,0,'SD',20000,'D:/temp/xxx.mat',200);
    res5_05 = detection2(st3,p);%include non user interfaced option
    
    save([dir(i,:),'cy5.mat'], 'res5_3', 'res5_05');
    clear('st3','res5_05','p');    
    
    %% colocalize the spots in channels 3 and 5  
    fp3sorted = sort_positions_arrays_to_match_spots(res5_3.final_pix{1},res3_3.final_pix{1});
    fp3sorted_clean = fp3sorted(1:size(res5_3.final_pix{1},1),:);
    clear('fp3sorted');
    n3matched = nnz(fp3sorted_clean(:,4));
    medianI3matched = median(nonzeros(fp3sorted_clean(:,4)));
    SDI3matched = std(nonzeros(fp3sorted_clean(:,4)));
    
    %make and save figure cy5 vs cy3 signal per particle 
    make_figure_mRNA_signal_Cy5_vs_Cy3(res5_3,res3_3,fp3sorted_clean,dir(i,:));
    
    %colocalize Cy5 and Cy3 signal to Cy3pt5 to get the transcription sites        
    spots3 = sort_positions_arrays_to_match_spots(res35_6.final_pix{1},res3_3.final_pix{1});
    spots5 = sort_positions_arrays_to_match_spots(res35_6.final_pix{1},res5_3.final_pix{1});
    spots35 = res35_6.final_pix{1};
    spots3(:,7) = i; %number of image in 7th column
    spots5(:,7) = i;
    spots35(:,7) = i;
    
    %save all results in a file
    medianI3 = median(res3_3.final_pix{1}(:,4));
    SDI3 = std(res3_3.final_pix{1}(:,4));
    n3 = size(res3_3.final_pix{1},1);
    medianI5 = median(res5_3.final_pix{1}(:,4));
    SDI5 = std(res5_3.final_pix{1}(:,4));
    n5 = size(res5_3.final_pix{1},1);
    save([dir(i,:),'allspots.mat'], 'spots3', 'spots5', 'spots35',...
        'medianI3','SDI3','n3','medianI5','SDI5','n5',...
        'medianI3matched','SDI3matched','n3matched');
    
    clear('spots3', 'spots5', 'spots35',...
        'medianI3','SDI3','n3','medianI5','SDI5','n5',...
        'medianI3matched','SDI3matched','n3matched',...
        'res35_6','res3_3','res5_3','fp3sorted_clean');
end
clear('dir','nst1','nst2','nst3');

end



