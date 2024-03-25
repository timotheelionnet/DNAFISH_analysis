%compile_2_color_analysis

%% load loc files for the 3 datasets (each dataset 2 color)

c1files = {'/Users/lionnett/Documents/Salina Long/two_color/3color_tim_1to1_cropped_datasets/C1-comp_tim2c_1to1_r1s3_lh_all_recon_trim.loc3';...
    '/Users/lionnett/Documents/Salina Long/two_color/3color_tim_1to1_cropped_datasets/C1-comp_tim2c_1to1_r1s4_rh_opticlobe_recon_trim_pos0-1.loc3';...
    '/Users/lionnett/Documents/Salina Long/two_color/3color_tim_1to1_cropped_datasets/C1-comp_tim2c_1to1_r1s5_lh_pdfcell_recon_trim-1.loc3'...
    };

%generate channel 2 dataset names by replacing the C1 string by C2
c2files = {'';'';'';};
for i=1:numel(c1files)
    c2files{i} = strrep(c1files{i},'C1','C2');
    
end

%load all the loc files
for i=1:numel(c1files)
    loc1{i} = dlmread(c1files{i});
    loc2{i} = dlmread(c2files{i});
end

%% parameters
%bin size for intensity histograms
binsize1 = 2000;
binsize2 = 2000;

%maximum displayed x cvalue for histograms
xmax1_for_plot = 5e6;
xmax2_for_plot = 2e6;

%guess fit values for 2 component lognormal function 
% y = A1 * exp( -(log(x) - I01).^2 / (2*s01^2) ) +  A2 * exp( -(log(x) - I02).^2 / (2*s02^2) )

%these guesses and fit bounds inputs correspond to the values read on the axis (will be log-
%transformed later)

%channel 1
A1guess1 = 350;
I01guess1 = 0.7e5;
s01guess1 = 1.5e4;

A2guess1 = 50;
I02guess1 = 1.7e5;
s02guess1 = 5e5;

Coeffs1 =   [A1guess1,  A2guess1,   I01guess1,  I02guess1,  s01guess1,  s02guess1];
lb1 =       [150,       5,          5e4,        1e5,        0,          0];
ub1 =       [1000,      70,         1e5,        4e5,        5e4,        5e5];

%channel 2
A1guess2 = 350;
I01guess2 = 1.4e5;
s01guess2 = 2e4;

A2guess2 = 50;
I02guess2 = 2.7e5;
s02guess2 = 5e5;

Coeffs2 =   [A1guess2,  A2guess2,   I01guess2,  I02guess2,  s01guess2,  s02guess2];
lb2 =       [150,       5,          1e5,        1.5e5,      0,          0];
ub2 =       [700,       50,         1.5e5,      5e5,        5e4,        5e6];

%log-transform guesstimates channel 1
lb1(3:end) = log(lb1(3:end));
ub1(3:end) = log(ub1(3:end));
Coeffs1(3:end) = log(Coeffs1(3:end));

%log-transform guesstimates channel 2
lb2(3:end) = log(lb2(3:end));
ub2(3:end) = log(ub2(3:end));
Coeffs2(3:end) = log(Coeffs2(3:end));  
%% plot Intensity histograms
%initialize figure for channel 1 datasets
fh1 = figure('Name','Channel 1 data'); 
height_offset = 0.1;
plot_height = (1 - 4*height_offset)/3;

axes1{1} = axes('Parent',fh1,'Position',[0.2,3*height_offset+2*plot_height,0.6,plot_height],'FontSize',12,'XScale','log');
hold(axes1{1},'on');


axes1{2} = axes('Parent',fh1,'Position',[0.2,2*height_offset+plot_height,0.6,plot_height],'FontSize',12,'XScale','log');
hold(axes1{2},'on');


axes1{3} = axes('Parent',fh1,'Position',[0.2,height_offset,0.6,plot_height],'FontSize',12,'XScale','log');
hold(axes1{3},'on');

%initialize figure for channel 2 datasets
fh2 = figure('Name','Channel 2 data'); 
height_offset = 0.1;
plot_height = (1 - 4*height_offset)/3;

axes2{1} = axes('Parent',fh2,'Position',[0.2,3*height_offset+2*plot_height,0.6,plot_height],'FontSize',12,'XScale','log');
hold(axes2{1},'on');


axes2{2} = axes('Parent',fh2,'Position',[0.2,2*height_offset+plot_height,0.6,plot_height],'FontSize',12,'XScale','log');
hold(axes2{2},'on');


axes2{3} = axes('Parent',fh2,'Position',[0.2,height_offset,0.6,plot_height],'FontSize',12,'XScale','log');
hold(axes2{3},'on');


for i=1:numel(loc1)
    %generate histogram data channel 1
    [ydata1{i},xdata1{i}] = hist(loc1{i}(:,4),0:binsize1:max(loc1{i}(:,4)));
    
    %plot hist raw data channel 1
    bar(xdata1{i},ydata1{i},'DisplayName','data','Parent',axes1{i},'LineWidth',2,...
    'FaceColor',[0.6 0.6 0.6],...
    'EdgeColor',[0.6 0.6 0.6]);
    axes(axes1{i});
    xlim([5e4 xmax1_for_plot]);
    ylim([0 max(ydata1{i}(2:end))]);
    
    %generate histogram data channel 2
    [ydata2{i},xdata2{i}] = hist(loc2{i}(:,4),0:binsize2:max(loc2{i}(:,4)));
    
    %plot hist raw data channel 2
    bar(xdata2{i},ydata2{i},'DisplayName','data','Parent',axes2{i},'LineWidth',2,...
    'FaceColor',[0.6 0.6 0.6],...
    'EdgeColor',[0.6 0.6 0.6]);
    axes(axes2{i});
    xlim([5e4 xmax2_for_plot]);
    ylim([0 max(ydata2{i}(2:end))]);
end

%% fit all histograms to 2-components lognormal (and plot)
y1bg = [];
y1sig = [];
y1tot = [];

y2bg = [];
y2sig = [];
y2tot = [];

for i=1:numel(loc1)
    
    %fit channel 1 data to 2-component lognormal
    options = optimset('TolX',.001,'MaxIter',1000);
    [Coeffs1(1) idx] = max(ydata1{i}(2:end));
    Coeffs1(3) = log(xdata1{i}(idx+1));
    [Coeffsout1{i}, resnorm]=lsqcurvefit(@multigaussian2, Coeffs1, log(xdata1{i}),ydata1{i},lb1,ub1,options);

    y1bg{i} = Coeffsout1{i}(1)*exp( -(log(xdata1{i}) - Coeffsout1{i}(3)).^2./(2*Coeffsout1{i}(5)^2));
    y1sig{i} = Coeffsout1{i}(2)*exp( -(log(xdata1{i}) - Coeffsout1{i}(4)).^2./(2*Coeffsout1{i}(6)^2));
    y1tot{i} = y1bg{i}+y1sig{i};
    
    %fit channel 2 data to 2-component lognormal
    options = optimset('TolX',.001,'MaxIter',1000);
    [Coeffs2(1) idx] = max(ydata2{i}(2:end));
    Coeffs2(3) = log(xdata2{i}(idx));
    [Coeffsout2{i}, resnorm]=lsqcurvefit(@multigaussian2, Coeffs2, log(xdata2{i}),ydata2{i},lb2,ub2,options);

    y2bg{i} = Coeffsout2{i}(1)*exp( -(log(xdata2{i}) - Coeffsout2{i}(3)).^2./(2*Coeffsout2{i}(5)^2));
    y2sig{i} = Coeffsout2{i}(2)*exp( -(log(xdata2{i}) - Coeffsout2{i}(4)).^2./(2*Coeffsout2{i}(6)^2));
    y2tot{i} = y2bg{i}+y2sig{i};
    
    %plot results channel 1
    plot(axes1{i},xdata1{i},y1bg{i},'LineWidth',2);
    plot(axes1{i},xdata1{i},y1sig{i},'LineWidth',2);
    plot(axes1{i},xdata1{i},y1tot{i},'LineWidth',2);
    
    %plot results channel 2
    plot(axes2{i},xdata2{i},y2bg{i},'LineWidth',2);
    plot(axes2{i},xdata2{i},y2sig{i},'LineWidth',2);
    plot(axes2{i},xdata2{i},y2tot{i},'LineWidth',2);
end

%% compute Jaccard and sort out spots based on Jaccard for each channel
jacc1 = [];
jacc2 = [];
locsort1 = [];
locsort2 = [];
for i=1:numel(loc1)
    for j=1:numel(xdata1{i})
        jacc1{i}(j) = sum(y1sig{i}(j:end))/(sum(y1sig{i}(:)) + sum(y1bg{i}(j:end)));
    end
    
    for j=1:numel(xdata2{i})
        jacc2{i}(j) = sum(y2sig{i}(j:end))/(sum(y2sig{i}(:)) + sum(y2bg{i}(j:end)));
    end
    
    [maxjac1(i),thresh1(i)] = max(jacc1{i});
    thresh1(i) = xdata1{i}(thresh1(i));
    [maxjac2(i),thresh2(i)] = max(jacc2{i});
    thresh2(i) = xdata2{i}(thresh2(i));
    
    locsort1{i} = loc1{i}(loc1{i}(:,4)>thresh1(i),:);
    locsort2{i} = loc2{i}(loc2{i}(:,4)>thresh2(i),:);
end


%% compute offset in x,y,z
f_dist_offset = figure; 
hold;
offset_binsize = 20; %in nm

order    = 5;
fcutlow  = 0.00001;
fcuthigh = 0.2;
[b,a]    = butter(order,[fcutlow,fcuthigh]*2, 'bandpass');

xvox = 88;
yvox = 88;
zvox = 100;

for i=1:numel(locsort1)
    % compute distance matrices
    n1 = size(locsort1{i},1);
    n2 = size(locsort2{i},1);

    %between the two channels
    dx12{i} = xvox*(repmat(locsort1{i}(:,1),1,n2) - repmat(locsort2{i}(:,1)',n1,1));
    dy12{i} = yvox*(repmat(locsort1{i}(:,2),1,n2) - repmat(locsort2{i}(:,2)',n1,1));
    dz12{i} = zvox*(repmat(locsort1{i}(:,3),1,n2) - repmat(locsort2{i}(:,3)',n1,1));
    dr12{i} = sqrt(dx12{i}.^2 + dy12{i}.^2 + dz12{i}.^2);

    % plot absolute offset in each dimensions
    x_edges = -1000:offset_binsize:1000;
    [nx,x_edges] = histcounts(dx12{i}(:),x_edges);
    [ny,y_edges] = histcounts(dy12{i}(:),x_edges);
    [nz,z_edges] = histcounts(dz12{i}(:),x_edges);

    
    plot((x_edges(2:end) + x_edges(1:end-1))/2,smooth(nx/sum(nx),5),'DisplayName',['dx12, ', num2str(i)]);
    plot((y_edges(2:end) + y_edges(1:end-1))/2,smooth(ny/sum(ny),5),'DisplayName',['dy12, ',num2str(i)]);
    plot((z_edges(2:end) + z_edges(1:end-1))/2,smooth(nz/sum(nz),5),'DisplayName',['dz12, ',num2str(i)]);
end

%% plot pair distribution
figure;
hold;
r_edges = 0:offset_binsize:10000;
for i=1:numel(locsort1)
    [nr{i},r_edges] = histcounts(dr12{i}(:),r_edges);
    denom = pi*(r_edges(2:end).^2 - r_edges(1:end-1).^2);
    plot((r_edges(2:end) + r_edges(1:end-1))/2,nr{i}./(denom*sum(nr{i})),'DisplayName',['dr12, ',num2str(i)]);
end

%%
%based on visual inspection
dist_thresh = [230,230,290];

for i=1:numel(locsort1)
    n1 = size(locsort1{i},1);
    n2 = size(locsort2{i},1);
    
    Idx{i}=zeros(n1,3);
    for j=1:size(dr12{i},1)
        
        [Idx{i}(j,2),Idx{i}(j,1)] = min(dr12{i}(j,:));
        
        if(Idx{i}(j,2)>dist_thresh(i))
            Idx{i}(j,3) = 0;
        else
            if ~ismember(Idx{i}(j,1), Idx{i}(:,3))
                Idx{i}(j,3) = Idx{i}(j,1);
            end
        end
    end
end


