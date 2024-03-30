% cell mask maker for pipe

maxZ_dir = char("D:\Finn\rep_nuclei_images_cellpose\processed");
my_fl = getFileListFromDir(maxZ_dir, 0, 0);

save_dir = 'D:\Finn\rep_nuclei_images_cellpose\processed\masks';

avg_nuke_size = 60; % in pixels

cell_thresh = 0;

for i = 1 :numel(my_fl)

    cur_zmax_path = my_fl{i};

    mask = generate_cellpose_mask(cur_zmax_path, avg_nuke_size , cell_thresh, save_dir);

    if mod(i, 10) ~=0 % lets keep evey 10th figrue for insepction

        close
    end


end

