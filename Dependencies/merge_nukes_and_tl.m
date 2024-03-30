config_fname = "H:\Finn_DriveH\20240226_PPARG_Staggered_Timecourse\adipo_rep2\adipo_picking_rep2\adipo_picking_rep2.ini";
config_fname = char(config_fname);

myParams = readParamsFromIniFile(config_fname);

res_dir = myParams.Output.outFolder;

fs = fileSet;

% populate object properties from ini config file
fs.buildFileSetFromConfig(config_fname);
fs.toNumeric({'FOV'},1);

%% loop thru and overlay

fov_vals = unique(fs.fList.FOV);

for i = 1:numel(fov_vals)
% for i =1:1

    cur_fov = fov_vals(i);

    cur_nuke_path = string(fs.getFileName({'FOV'}, {cur_fov}, {'nukes'}));

    cur_tl_path = string(fs.getFileName({'FOV'}, {cur_fov}, {'tl'}));

    cur_nuke = imread(cur_nuke_path);

    cur_tl = imread(cur_tl_path);

    cur_tl = cur_tl;
    
    figure
    comp = imfuse(cur_nuke, cur_tl, 'falsecolor','Scaling','independent','ColorChannels',[1 2 0]);

    imshow(comp)

    f_name = fullfile(res_dir, 'tl_nuke_comp_' + string(i) +'.png');



    imwrite(comp, f_name)





end
