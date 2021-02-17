
%inputs
%res is the result from AIRLOCALIZE
%a threshold for signal
signal_threshold = 90000;

%size of stack in ImageJ order
stack_size = [414,412,427];

%input tif file location (ensure 1-channel image)
tifdata = '/Users/lionnett/Documents/Salina Long/two_color/101716/C1-comp_tim2c_2to1_r1s6_rh_all_recon_trim.tif';

%histogram binsize
binsize = 2000;
%% histogram data of spots Intensities including background detections
[ydata,xdata] = hist(res.final_pix(:,4),0:binsize:max(res.final_pix(:,4)));
figure;
bar(xdata,ydata,'DisplayName','data','LineWidth',2,...
    'FaceColor',[0.6 0.6 0.6],...
    'EdgeColor',[0.6 0.6 0.6]);%raw data
%% s
spots_artefacts = res.final_pix(res.final_pix(:,4)<0,:);
spots_artefacts(:,4) = 1000;

spots_bg = res.final_pix(res.final_pix(:,4)>0 & res.final_pix(:,4)<signal_threshold,:);
spots_sig = res.final_pix(res.final_pix(:,4)>signal_threshold,:);

stack_sizeMatlab = [fliplr(stack_size(1:2)),stack_size(3)];

stack_artefacts = generate_gauss_spots_in_stack(spots_artefacts,[1,1,1],stack_sizeMatlab);
stack_bg = generate_gauss_spots_in_stack(spots_bg,[1,1,1],stack_sizeMatlab);
stack_sig = generate_gauss_spots_in_stack(spots_sig,[1,1,1],stack_sizeMatlab);

%stack_artefacts = permute(stack_artefacts,[2,1,3]);
%stack_bg = permute(stack_bg,[2,1,3]);
%stack_sig = permute(stack_sig,[2,1,3]);

[f,n,e] = fileparts(tifdata);

save_as_tiff(stack_artefacts, fullfile(f,[n,'_artefact_spots.tif']));
save_as_tiff(stack_bg, fullfile(f,[n,'_bg_spots.tif']));
save_as_tiff(stack_sig, fullfile(f,[n,'_sig_spots.tif']));

%%