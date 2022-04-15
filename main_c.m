set(0,'defaultfigurecolor',[1 1 1])
clear, close all
addpath Sub
%rng('default')

%% DATA COLLECTION --------------------------------------------------------
% Topology parameters
K_class = 2;
sample_rate = 1;
meter_pixel = 9;
map_height = 50;    % Drone height for air-to-ground power map
sample_height = 50; % Drone height for learning (Learning stage)
neighbour_mat = [0 0; -1 0; 0 -1; 1 0; 0 1] * 39; %
neighbour_weight = [0.6; 0.1; 0.1; 0.1; 0.1]; % Sum-to-one
test_positions = round(1:100);
noise_guass = 3; % dB std
n_ue = 50; % [8*8 9*10 11*11 13*13 15*16 18*18 21*22 25*26 30*30 35*36]*8
n_uav = 50 * length(sample_height);    % 15:50:9 50:200:5 100:400
n_mea = n_ue * n_uav;
n_uav = round(n_uav / length(sample_height));
residual_on = false;
residsec_on = false;
residthi_on = false;
residfou_on = false;
residfiv_on = true;
segment_on = false;
segresid_on = false;
plot_on = false;
kmeans_on = false;

% DATA = load('radiomap_shanghai100tx_2.5GHz.mat');
% DATA = load('radiomap_simulated100tx_noise.mat');
% DATA = load('radiomap_shanghai115tx_28GHz.mat');
% DATA = load('radiomap_simulated100tx_nonoise');
DATA = load('radiomap_simulated100tx_3class.mat');
pos_ue = DATA.PosUE;
Xdata = zeros(n_mea, 6); 
Gdata = zeros(n_mea, 1);
Ddata = zeros(n_mea, 1);
Gtrue = cell(n_ue, 1);
cnt = 0;

figure(1)
title("Users' locations")
hold on
pos_ue_all = zeros(size(pos_ue));
userids = sort(randperm(length(pos_ue), n_ue));
for id = 1:length(userids)
    % i must increase
    i = userids(id);
    % Training data
    for j = 1:length(sample_height)
        uav_height = sample_height(j);
        x = pos_ue(i,1); y = pos_ue(i,2); z = pos_ue(i,3);
        I = (DATA.RadioMap(:, 6) == uav_height) ...
            & (DATA.RadioMap(:, 1) == x) ...
            & (DATA.RadioMap(:, 2) == y)...
            & (DATA.RadioMap(:, 3) == z);
        Rm2D = DATA.RadioMap(I, :);
        Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
        [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
        Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
        position = [x-min(Xvec)+1 y-min(Yvec)+1 z];
        Xmat = Xmat - min(Xvec) + 1;
        Ymat = Ymat - min(Yvec) + 1;
        % Data pre-processing
        Zmat(isnan(Zmat)) = min(min(Zmat));
        Zmat = Zmat + randn(size(Zmat)) * noise_guass;
        Xvec = Xmat(:);
        Yvec = Ymat(:);
        Zvec = Zmat(:);
        dist = log10(vecnorm([Xvec Yvec uav_height*ones(length(Xvec),1)] - ...
            ones(length(Xvec),3).*position, 2, 2));
        position = [ceil(position(1:2)) z];
        if ~exist('idx', 'var')
            %idx = 1:ceil(length(Xvec)/nDrone):length(Xvec);
            idx = randperm(length(Xvec), n_uav);
        end
        uav_positon = ceil([Xvec(idx) Yvec(idx)]);
        uav_num = length(uav_positon);
        cnt = cnt + uav_num;
        Xdata(cnt-uav_num+1:cnt, :) = [uav_positon ...
            uav_height.*ones(uav_num,1) ones(uav_num,3).*position];
        Gdata(cnt-uav_num+1:cnt, :) = Zvec(idx);
        Ddata(cnt-uav_num+1:cnt, :) = dist(idx);

        figure(1)
        plot(position(1), position(2), '.');
        text(position(1), position(2), string(i));

        pos_ue_all(i, :) = position;
    end
    % Test data
    uav_height = map_height;
    x = pos_ue(i,1); y = pos_ue(i,2); z = pos_ue(i,3);
    I = (DATA.RadioMap(:, 6) == uav_height) ...
        & (DATA.RadioMap(:, 1) == x) ...
        & (DATA.RadioMap(:, 2) == y)...
        & (DATA.RadioMap(:, 3) == z);
    Rm2D = DATA.RadioMap(I, :);
    Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
    [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
    Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
    position = [x-min(Xvec)+1 y-min(Yvec)+1 z];
    Xmat = Xmat - min(Xvec) + 1;
    Ymat = Ymat - min(Yvec) + 1;
    % Data pre-processing
    Xvec = Xmat(:);
    Yvec = Ymat(:);
    Zvec = Zmat(:);
    Zmat(isnan(Zmat)) = min(Zvec);
    Gtrue{i} = Zmat;
end
lenx = ceil(max(Xvec)/meter_pixel);
leny = ceil(max(Yvec)/meter_pixel);
n_mea = cnt;
Xdata = Xdata(1:cnt, :); 
Ydata = Gdata(1:cnt, :);
Ddata = Ddata(1:cnt, :);
maps.BldMapZ = zeros(lenx, leny);
maps.BldPosMat = ones(lenx, leny);
maps.FolPosMat = zeros(lenx, leny);
maps.meterPerPixel = meter_pixel;

% Generate collinear obstacles set
obstacles = (1:lenx*leny)';
nObst = length(obstacles);
cols = cell(nObst, 2);   % collinear measurements of each obstacle
covB = cell(n_mea, 2);
for id = 1:n_mea
    uav_pos_meter = Xdata(id, 1:3);
    uav_pos_pixel = [floor(uav_pos_meter(1:2)/meter_pixel)+1, uav_pos_meter(3)];
    ue_pos_meter = Xdata(id, 4:6);
    ue_pos_pixel = [floor(ue_pos_meter(1:2)/meter_pixel)+1, ue_pos_meter(3)];
    [cov_bld, cov_z] = covBldZ(uav_pos_pixel, ue_pos_pixel, lenx);
    cov_bld = cov_bld(cov_z > 0 & cov_z <= map_height);
    cov_z = cov_z(cov_z > 0 & cov_z <= map_height);
    covB{id, 1} = cov_bld;
    covB{id, 2} = cov_z;
    for i = 1:length(cov_bld)
        j = cov_bld(i);
        ib = find(obstacles == j, 1, 'first');
        if ~isempty(ib)
            if ~isempty(cols{ib, 1})
                cols{ib, 1} = [cols{ib, 1} id];
                cols{ib, 2} = [cols{ib, 2} cov_z(i)];
            else
                cols{ib, 1} = id;
                cols{ib, 2} = cov_z(i);
            end
        end
    end
end

%% LOCAL POLY RECONSTRUCTION ----------------------------------------------
MAXLOOP = 9;
tolerance = 1e-3;
metrics = cell(MAXLOOP, 6);
maps.droneHeightMap = map_height;
maps.neighbourMat = neighbour_mat;
maps.neighbourWeight = neighbour_weight;

% Initialization
R.X = Xdata;
P = polyfit(Ddata, Ydata, 1);
delta = Ydata - (P(1) * Ddata + P(2));
w = delta > 0;
W = [w w ~w ~w];
A = [Ddata ones(n_mea, 1)];
A_ = W .* [A A];
X = (A_' * A_) \ A_' * Ydata;
A = repmat(A, 1, K_class+1);
R.Z1 = zeros(n_mea, K_class+1);
R.Alpha = zeros(1, K_class+1);
R.Beta = zeros(1, K_class+1);
for k = 1:K_class+1
    R.Alpha(k) = X(1) - (X(1) - X(3)) * ((k-1)/K_class);
    R.Beta(k) = X(2) - (X(2) - X(4)) * ((k-1)/K_class);
    R.Z1(:, k) = -abs(Ydata-(R.Alpha(k)*Ddata+R.Beta(k)));
end

% -- coarse resolution --
dsfactor = 4;       % Downsampling factor
[small_maps, S, Hs] = downsampleMaps(maps, dsfactor, R);
R.Hs = Hs; R.S = S;
smap_hat = ones([K_class numel(small_maps.BldMapZ)]) * map_height;

for i = 1:MAXLOOP
    [smap_hat, W] = optimizeH(R, small_maps, smap_hat);
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    R.Alpha = X(1:2:end);
    R.Beta = X(2:2:end);
    for k = 1:K_class+1
        R.Z1(:, k) = -abs(Ydata-(R.Alpha(k)*Ddata+R.Beta(k)));
    end
end

% -- medium resolution --
dsfactor2 = 2;  
[medium_maps, S, Hs] = downsampleMaps(maps, dsfactor2, R);
R.Hs = Hs; R.S = S;
mmap_hat = zeros([K_class numel(medium_maps.BldMapZ)]);
for k = 1:K_class
    upsample_matrix = repelem(reshape(smap_hat(k, :), ...
        size(small_maps.BldMapZ)), round(dsfactor/dsfactor2), ...
        round(dsfactor/dsfactor2));
    mmap_hat(k, :) = vec(upsample_matrix(1:size(medium_maps.BldMapZ, 1),...
        1:size(medium_maps.BldMapZ, 2)));
end

for i = 1:MAXLOOP
    [mmap_hat, W] = optimizeH(R, medium_maps, mmap_hat);
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    R.Alpha = X(1:2:end);
    R.Beta = X(2:2:end);
    for k = 1:K_class+1
        R.Z1(:, k) = -abs(Ydata-(R.Alpha(k)*Ddata+R.Beta(k)));
    end
end

% -- fine resolution (3 meter) --
dsfactor3 = 1;
[maps, S, Hs] = downsampleMaps(maps, dsfactor3, R);
R.Hs = Hs; R.S = S;
map_hat = zeros([K_class numel(maps.BldMapZ)]);
for k = 1:K_class
    upsample_matrix = repelem(reshape(mmap_hat(k, :), ...
        size(medium_maps.BldMapZ)), round(dsfactor2/dsfactor3), ...
        round(dsfactor2/dsfactor3));
    map_hat(k, :) = vec(upsample_matrix(1:size(maps.BldMapZ, 1),...
        1:size(maps.BldMapZ, 2)));
end

h  = waitbar(0, 'Estimate building heights');
for i = 1:MAXLOOP
    waitbar(i / MAXLOOP, h);
    [map_hat, W] = optimizeH(R, maps, map_hat);
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    R.Alpha = X(1:2:end);
    R.Beta = X(2:2:end);
    for k = 1:K_class+1
        R.Z1(:, k) = -abs(Ydata-(R.Alpha(k)*Ddata+R.Beta(k)));
    end
end
close(h);

% hf = showmap(BldMapHat + FolMapHat, Maps.meterPerPixel, 8);
% title('Esitmated building and foliage height');

if kmeans_on == true
    c_i = D1 ./ abs((X(1)*Ddata+X(2)) - (X(3)*Ddata+X(4)));
    obs_clu = zeros(nObst, map_height);
    for i = 1:map_height
        for j = 1:nObst
            obs_clu(j, i) = mean(c_i(cols{j, 1}(ceil(cols{j, 2}) == i)));
        end
    end
    obs_clu(isnan(obs_clu)) = 0;
    bld_ind = zeros(lenx, leny);
    bld_ind(maps.BldPosMat > 0) = kmeans(obs_clu, 5);
end

figure, hold on
for k = 1:K_class+1
    plot(Ddata(W(:, 2*k) > 1/(K_class+1)), ...
        Ydata(W(:, 2*k) > 1/(K_class+1)), '.', 'Color', rand([1 3]));
end
figure, hold on
for k = K_class+1:-1:1
    plot(Ddata(W(:, 2*k) > 1/(K_class+1)), ...
        Ydata(W(:, 2*k) > 1/(K_class+1)), '.', 'Color', rand([1 3]));
end
if isfield(DATA, 'BldPosMat')
    k = K_class;
    map_hat(k, DATA.BldPosMat(round(1:DATA.lenX/lenx:DATA.lenX), ...
        round(1:DATA.lenY/leny:DATA.lenY)) < 1) = 0;
end
heights = reshape(map_hat(K_class, :), size(maps.BldMapZ));
for i = 1:length(pos_ue_all)
    x = round(pos_ue_all(i, 1) / meter_pixel);
    y = round(pos_ue_all(i, 2) / meter_pixel);
    if x == 0 || y == 0, continue, end
    heights(x, y) = 0;
end
% heights(heights >= droneHeightMap | heights == droneHeightMap/2) = 0;
% indicator = heights >= UserHeight;
% heights = medfilt1(heights .* (heights <= DroneHeight)) .* indicator;
figure; showmap(heights, meter_pixel);
title(sprintf('Height from L model at %d meters', round(map_height)));

uav_num = uav_num * length(sample_height);

%% RADIOMAP AND PERFORMANCE -----------------------------------------------
% Channel model 
C.A1 = X(1); C.B1 = X(2); C.S1 = 0;
C.A2 = 0; C.B2 = 0; C.S2 = 0; 
C.A3 = X(3); C.B3 = X(4); C.S3 = 0;
channel_model.C = C;
channel_model.los_nlos_trans = 0;
channel_model.noise = 0;

metrics_mse = []; metrics_mae = [];
XYkr = [vec(meshgrid(1:lenx, 1:leny)'), vec(meshgrid(1:leny, 1:lenx)), ...
    ones(lenx * leny, 1) * map_height / meter_pixel];
Xann = [vec(meshgrid(1:lenx, 1:leny)'), vec(meshgrid(1:leny, 1:lenx)) ...
    * meter_pixel, ones(lenx * leny, 1) * map_height];
net = fitnet([32 32 16 16 8 8 4 4]); %fitnet([(4:8).^2 (8:-1:2).^2]);
net = train(net, [Xdata Ddata X(1)*Ddata+X(2) X(3)*Ddata+X(4)]', Ydata', ...
    'useParallel', 'no', 'showResources', 'no');

j = 1;
for i = 1:length(pos_ue)
    if pos_ue_all(i, 1) == 0 || pos_ue_all(i, 2) == 0, continue, end
    if ~ismember(i, test_positions), j = j + 1; continue, end
    Gtru = Gtrue{i}';
    % Radiomap reconstruction
    Xdata_i = [Xdata((j-1)*uav_num+1:j*uav_num,1:2)/meter_pixel, ...
        Xdata((j-1)*uav_num+1:j*uav_num,3)/meter_pixel];
    Ydata_i = Ydata((j-1)*uav_num+1:j*uav_num);
    idx = randperm(length(Xdata_i), round(length(Xdata_i) / sample_rate));
    Xdata_i = Xdata_i(idx, :); Ydata_i = Ydata_i(idx, :);
    d_variogram = variogram(Xdata_i(:, 1:2), Ydata_i);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    if n_mea > 1e9
        idx = randperm(length(Xdata_i), round(length(Xdata_i) / sample_rate));
        Gkri = kriging(vstruct,Xdata_i(idx,:),false,Ydata_i(idx),XYkr,false);
    else
        idx = randperm(length(Xdata), round(length(Xdata) / sample_rate));
        Gkri = kriging(vstruct, Xdata(idx,:)/meter_pixel, false, Ydata(idx), ...
        [XYkr ones(length(XYkr), 3).*pos_ue_all(i, :)/meter_pixel], false);
    end
    Xann_upos = ones(length(XYkr), 3).*pos_ue_all(i, :);
    Xann_dist = log10(vecnorm(Xann - Xann_upos, 2, 2));
    Xann_features = [Xann Xann_upos ...
        Xann_dist X(1)*Xann_dist+X(2) X(3)*Xann_dist+X(4)]';
    Gann = net(Xann_features);
    Gann = reshape(Gann, [lenx leny]);
    Gkri = reshape(Gkri, [lenx leny]);
    Gknn = powerMapReconKNN(R,Ydata,pos_ue_all(i, :),map_height,maps);
    Gbld = powerMapReconBld(R,heights,heights,pos_ue_all(i, :),map_height,maps);
    Gbdk = powerMapReconKBld(R, map_hat, pos_ue_all(i, :), map_height, maps);
    % Radiomap reconstruction with residual
    if residual_on == true
    [~, id] = min(abs(sample_height - map_height));
    droneHeightRes = sample_height(id)/meter_pixel;
    ResX = floor(Xdata_i(:, 1:2)) + 1;
    ResX = ResX(Xdata_i(:, 3) == droneHeightRes, :);
    ResG = Ydata_i(Xdata_i(:, 3) == droneHeightRes);
    Gbld = powerMapReconBld(R,heights,heights,pos_ue_all(i, :),map_height,maps);
    Gbld = imgaussfilt(Gbld, 1);
    for k = 1:length(ResX)
        ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
    end
    d_variogram = variogram(ResX, ResG);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    Rkri = kriging(vstruct, ResX(:,1), ResX(:,2), ...
        ResG, meshgrid(1:lenx, 1:leny)', meshgrid(1:leny, 1:lenx));
    Gbld = Gbld + reshape(Rkri, [lenx leny]);
    end
    if residsec_on == true
    ResX = Xdata / meter_pixel;
    P = (X(1)*Ddata + X(2) - Ydata)./(X(1)*Ddata + X(2) - X(3)*Ddata - X(4));
    P(P < 0) = 0; P(P > 1) = 1;
    ResG = Ydata - w1.*(X(1)*Ddata + X(2)) - w2.*(X(3)*Ddata + X(4));
    d_variogram = variogram(Xdata((j-1)*uav_num+1:j*uav_num,1:2)/meter_pixel, ...
        Ydata((j-1)*uav_num+1:j*uav_num));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    ResRec = kriging(vstruct, ResX, false, ResG, ...
    [XYkr ones(length(XYkr), 3).*pos_ue_all(i, :)/meter_pixel], false);
    ResRec = reshape(ResRec, [lenx leny]);
    Gbld = Gbld + ResRec;
    end
    if residthi_on == true
    Rind = true(length(ResX), 1);
    for k = 1:length(ResX)
        ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
        Rind(k) = Gind(ResX(k, 1), ResX(k, 2));
    end
    d_variogram = variogram(ResX(Rind,:), ResG(Rind));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    LOSSeg = kriging(vstruct, ResX(Rind,1), ResX(Rind,2), ...
        ResG(Rind), meshgrid(1:lenx, 1:leny)', meshgrid(1:leny, 1:lenx));
    d_variogram = variogram(ResX(~Rind,:), ResG(~Rind));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    NLOSSeg = kriging(vstruct, ResX(~Rind,1), ResX(~Rind,2), ...
        ResG(~Rind), meshgrid(1:lenx, 1:leny)', meshgrid(1:leny, 1:lenx));
    Gbld = Gbld + LOSSeg .* Gind + NLOSSeg .* ~Gind;
    end
    if residfou_on == true
    ResY = []; ResR = [];
    for id = 1:length(sample_height)
        h = sample_height(id)/meter_pixel;
        ResX = [floor(Xdata_i(:, 1:2)) + 1 Xdata_i(:, 3)];
        ResX = ResX(Xdata_i(:, 3) == h, :);
        ResG = Ydata_i(Xdata_i(:, 3) == h);
        Gbld = powerMapReconBld(R,heights,heights,pos_ue_all(i, :),h,maps);
        Gbld = imgaussfilt(Gbld, 3);
        for k = 1:length(ResX)
            ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
        end
        x = ResX(:, 1) - pos_ue_all(i, 1)/meter_pixel;
        y = ResX(:, 2) - pos_ue_all(i, 2)/meter_pixel;
        h = h/meter_pixel - pos_ue_all(i, 3)/meter_pixel;
        l = vecnorm([x y], 2, 2);
        theta1 = atan2(y,x); theta2 = atan2(h,l);
        ResR = [ResR; min(theta1+pi,pi-theta1)/pi*180 ...
            min(theta2+pi,pi-theta2)/pi*180 sqrt(l.^2 + h^2)];
        ResY = [ResY; ResG];
    end
    d_variogram = variogram(ResR, ResY);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    Rkri = krigingR(vstruct, Xdata((j-1)*uav_num+1:j*uav_num,1:3)/ ...
        meter_pixel, ResY, XYkr, pos_ue_all(i, :)/meter_pixel);
    Gbld = powerMapReconBld(R,heights,heights,pos_ue_all(i, :),map_height,maps);
    Gbld = imgaussfilt(Gbld, 3) + reshape(Rkri, [lenx leny]);
    end
    if residfiv_on == true
    [~, id] = min(abs(sample_height - map_height));
    droneHeightRes = sample_height(id)/meter_pixel;
    ResX = floor(Xdata_i(:, 1:2)) + 1;
    ResX = ResX(Xdata_i(:, 3) == droneHeightRes, :);
    ResG = Ydata_i(Xdata_i(:, 3) == droneHeightRes);
    Gind = powMapRecBldSof(R,heights,[],pos_ue_all(i, :),map_height,maps,[]);
    Gind = imgaussfilt(Gind, 3);
    Glos = powMapRecBldSof(R,heights-realmax,[],pos_ue_all(i, :),map_height,maps);
    Gnlo = powMapRecBldSof(R,heights+realmax,[],pos_ue_all(i, :),map_height,maps);
    Gbld = Gind .* Glos + (1 - Gind) .* Gnlo;
    Gbdn = Gbld;
    for k = 1:length(ResX)
        ResG(k) = ResG(k) - Gbld(min(lenx, ResX(k, 1)), min(leny, ResX(k, 2)));
    end
    d_variogram = variogram(ResX, ResG);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    Rkri = kriging(vstruct, ResX(:,1), ResX(:,2), ...
        ResG, meshgrid(1:lenx, 1:leny)', meshgrid(1:leny, 1:lenx));
    Gbld = Gbld + reshape(Rkri, [lenx leny]);
    end
    % Radiomap reconstruction with segmentation
    if segment_on == true
    d_variogram = variogram(Xdata(w1, :)/meter_pixel, Ydata(w1));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    LOSSeg = kriging(vstruct,Xdata(w1, :)/meter_pixel, false, Ydata(w1), ...
        [XYkr ones(length(XYkr), 3).*pos_ue_all(i, :)/meter_pixel], false);
    d_variogram = variogram(Xdata(w2, :)/meter_pixel, Ydata(w2));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    NLOSSeg = kriging(vstruct,Xdata(w2, :)/meter_pixel, false, Ydata(w2), ...
        [XYkr ones(length(XYkr), 3).*pos_ue_all(i, :)/meter_pixel], false);
    LOSSeg = reshape(LOSSeg, [lenx leny]);
    NLOSSeg = reshape(NLOSSeg, [lenx leny]);
    Gind = powerMapReconBld(R,heights,heights,pos_ue_all(i, :),map_height,maps,1);
    Gbld = LOSSeg .* Gind + NLOSSeg .* ~Gind;
    end
    % Radiomap reconstruction with residual and segmentation
    if segresid_on == true
    Xlos = Xdata(w1, :) / meter_pixel; 
    Ylos = Ydata(w1) - Ddata(w1) * X(1) - X(2);
    Xnlos = Xdata(w2, :) / meter_pixel;
    Ynlos = Ydata(w2) - Ddata(w2) * X(3) - X(4);
    d_variogram = variogram(Xlos, Ylos);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    LOSSeg = kriging(vstruct, Xlos, false, Ylos, ...
        [XYkr ones(length(XYkr), 3).*pos_ue_all(i, :)/meter_pixel], false);
    d_variogram = variogram(Xnlos, Ynlos);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    NLOSSeg = kriging(vstruct, Xnlos, false, Ynlos, ...
        [XYkr ones(length(XYkr), 3).*pos_ue_all(i, :)/meter_pixel], false);
    LOSSeg = reshape(LOSSeg, [lenx leny]);
    NLOSSeg = reshape(NLOSSeg, [lenx leny]);
    Gind = powerMapReconBld(R,heights,heights,pos_ue_all(i, :),map_height,maps,1);
    Gbld = Gbld + LOSSeg .* Gind + NLOSSeg .* ~Gind;
    end
    % Performance and figures
    xm = length(Zmat);
    if xm <= lenx
        Gbld = Gbld(round(1:lenx/xm:lenx),round(1:lenx/xm:lenx));
        Gkri = Gkri(round(1:lenx/xm:lenx),round(1:lenx/xm:lenx));
        Gknn = Gknn(round(1:lenx/xm:lenx),round(1:lenx/xm:lenx));
        Gann = Gann(round(1:lenx/xm:lenx),round(1:lenx/xm:lenx));
        Gbdk = Gbdk(round(1:lenx/xm:lenx),round(1:lenx/xm:lenx));
    else
        Gtru = Gtru(round(1:xm/lenx:xm),round(1:xm/lenx:xm));
    end
    if plot_on == true
    figure(300 + i);
    subplot(2, 3, 1); showmap(Gkri, meter_pixel);
    subplot(2, 3, 2); showmap(Gknn, meter_pixel);
    subplot(2, 3, 3); showmap(Gann, meter_pixel);
    subplot(2, 3, 4); showmap(Gbdk, meter_pixel);
    subplot(2, 3, 5); showmap(Gbld, meter_pixel);
    subplot(2, 3, 6); showmap(Gtru, meter_pixel);
    % figure(400 + i); showmap(reshape(Rkri, [lenX lenY]), 9);
    disp_str = strcat('Position ',string(i),': ',string(mean(abs(vec(...
        Gtru-Gkri)))),' vs. ',string(mean(abs(vec(Gtru-Gknn)))),' vs. ',...
        string(mean(abs(vec(Gtru-Gbld)))),'; ',string(mse(Gtru-Gkri)),...
        ' vs. ',string(mse(Gtru-Gknn)),' vs. ',string(mse(Gtru-Gbld)));
    disp(disp_str); input('');
    end
    metrics_mse = [metrics_mse; mse(Gtru-Gkri) mse(Gtru-Gknn) ...
        mse(Gtru-Gann) mse(Gtru-Gbdk) mse(Gtru-Gbld)];
    metrics_mae = [metrics_mae; mean(abs(vec(Gtru-Gkri))) ...
        mean(abs(vec(Gtru-Gknn))) mean(abs(vec(Gtru-Gann))) ...
        mean(abs(vec(Gtru-Gbdk))) mean(abs(vec(Gtru-Gbld)))];
    j = j + 1;
end
disp([mean(metrics_mae, 1) mean(metrics_mse, 1) n_ue]);
