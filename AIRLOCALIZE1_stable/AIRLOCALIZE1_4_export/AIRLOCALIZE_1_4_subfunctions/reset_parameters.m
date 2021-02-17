function params = reset_parameters(params)
    % removing info specific to one image so that it does not affect the
    % treatment of the next image
    
    % ROI
    if isfield(params,'select_ROI')
        params = rmfield(params,'select_ROI');
    end
    if isfield(params,'Xmin')
        params = rmfield(params,'Xmin');
    end
    if isfield(params,'Ymin')
        params = rmfield(params,'Ymin');
    end
    if isfield(params,'Xmax')
        params = rmfield(params,'Xmax');
    end
    if isfield(params,'Ymax')
        params = rmfield(params,'Ymax');
    end
    
    %img size
    if isfield(params,'nx')
        params = rmfield(params,'nx');
    end
    if isfield(params,'ny')
        params = rmfield(params,'ny');
    end
    if isfield(params,'nz')
        params = rmfield(params,'nz');
    end
    
    %data
    if isfield(params,'data');
        params = rmfield(params,'data');
    end
    if isfield(params,'smooth');
        params = rmfield(params,'smooth');
    end
    
    %%current file name
    if isfield(params,'curfname')
        params = rmfield(params,'curfname');
    end
end