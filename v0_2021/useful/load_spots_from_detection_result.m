function spots = load_spots_from_detection_result(spots_filename)

[~,~,ext] = fileparts(spots_filename);
switch ext
    case '.det'
            res = load(spots_filename,'-mat');
            res = res.detection_result;
            spots = res.final_pix;   
    case '.loc3'
            spots = load(spots_filename);
    case '.loc'
            spots = load(spots_filename);
end

end