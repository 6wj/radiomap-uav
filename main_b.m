%%
for testPhaseNum = 1:1
set(0,'defaultfigurecolor',[1 1 1])
clear, close all
addpath Sub
%rng('default')

%% DATA COLLECTION --------------------------------------------------------
% Topology parameters
sampleRate = 1;
meterPerPixel = 9;
droneHeightMap = 50;    % Drone height for air-to-ground power map
droneHeightSample = 50:10:120; % Drone height for learning (Learning stage)
neighbourMat = [0 0]; %[0 0; -1 0; 0 -1; 1 0; 0 1] * 39; %
testPositions = round(1:100);
noiseGuass = 3; % dB std
nUser = 50; % [8*8 9*10 11*11 13*13 15*16 18*18 21*22 25*26 30*30 35*36]*8
nDrone = 50 * length(droneHeightSample);    % 15:50:9 50:200:5 100:400
nMeas = nUser * nDrone;
nDrone = round(nDrone / length(droneHeightSample));
residualOn = false;
residsecOn = false;
residthiOn = false;
residfouOn = false;
residfivOn = true;
segmentOn = false;
segresidOn = false;
plotOn = false;
kmeansOn = false;

DATA = load('radiomap_shanghai100tx_2.5GHz.mat');
% DATA = load('radiomap_shanghai115tx_28GHz.mat');
PosUE = DATA.PosUE;
Xdata = zeros(nMeas, 6); 
Gdata = zeros(nMeas, 1);
Ddata = zeros(nMeas, 1);
Gtrue = cell(nUser, 1);
cnt = 0;

figure(1)
title("Users' locations")
hold on
userPosAll = zeros(size(PosUE));
userids = sort(randperm(length(PosUE), nUser));
for id = 1:length(userids)
    % i must increase
    i = userids(id);
    % Training data
    for j = 1:length(droneHeightSample)
        DroneHeight = droneHeightSample(j);
        x = PosUE(i,1); y = PosUE(i,2); z = PosUE(i,3);
        I = (DATA.RadioMap(:, 6) == DroneHeight) ...
            & (DATA.RadioMap(:, 1) == x) ...
            & (DATA.RadioMap(:, 2) == y)...
            & (DATA.RadioMap(:, 3) == z);
        Rm2D = DATA.RadioMap(I, :);
        Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
        [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
        Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
        UserPos = [x-min(Xvec)+1 y-min(Yvec)+1 z];
        Xmat = Xmat - min(Xvec) + 1;
        Ymat = Ymat - min(Yvec) + 1;
        % Data pre-processing
        Zmat(isnan(Zmat)) = min(min(Zmat));
        Zmat = Zmat + randn(size(Zmat)) * noiseGuass;
        Xvec = Xmat(:);
        Yvec = Ymat(:);
        Zvec = Zmat(:);
        dist = log10(vecnorm([Xvec Yvec DroneHeight*ones(length(Xvec),1)] - ...
            ones(length(Xvec),3).*UserPos, 2, 2));
        UserPos = [ceil(UserPos(1:2)) z];
        if ~exist('idx', 'var')
            %idx = 1:ceil(length(Xvec)/nDrone):length(Xvec);
            idx = randperm(length(Xvec), nDrone);
        end
        DronePositon = ceil([Xvec(idx) Yvec(idx)]);
        DroneNum = length(DronePositon);
        cnt = cnt + DroneNum;
        Xdata(cnt-DroneNum+1:cnt, :) = [DronePositon ...
            DroneHeight.*ones(DroneNum,1) ones(DroneNum,3).*UserPos];
        Gdata(cnt-DroneNum+1:cnt, :) = Zvec(idx);
        Ddata(cnt-DroneNum+1:cnt, :) = dist(idx);

        figure(1)
        plot(UserPos(1), UserPos(2), '.');
        text(UserPos(1), UserPos(2), string(i));

        userPosAll(i, :) = UserPos;
    end
    % Test data
    DroneHeight = droneHeightMap;
    x = PosUE(i,1); y = PosUE(i,2); z = PosUE(i,3);
    I = (DATA.RadioMap(:, 6) == DroneHeight) ...
        & (DATA.RadioMap(:, 1) == x) ...
        & (DATA.RadioMap(:, 2) == y)...
        & (DATA.RadioMap(:, 3) == z);
    Rm2D = DATA.RadioMap(I, :);
    Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
    [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
    Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
    UserPos = [x-min(Xvec)+1 y-min(Yvec)+1 z];
    Xmat = Xmat - min(Xvec) + 1;
    Ymat = Ymat - min(Yvec) + 1;
    % Data pre-processing
    Xvec = Xmat(:);
    Yvec = Ymat(:);
    Zvec = Zmat(:);
    Zmat(isnan(Zmat)) = min(Zvec);
    Gtrue{i} = Zmat;
end
lenX = ceil(max(Xvec)/meterPerPixel);
lenY = ceil(max(Yvec)/meterPerPixel);
nMeas = cnt;
Xdata = Xdata(1:cnt, :); 
Ydata = Gdata(1:cnt, :);
Ddata = Ddata(1:cnt, :);
Maps.BldMapZ = zeros(lenX, lenY);
Maps.BldPosMat = ones(lenX, lenY);
Maps.FolPosMat = zeros(lenX, lenY);
Maps.meterPerPixel = meterPerPixel;

R.X = Xdata;
P = polyfit(Ddata, Ydata, 1);
BldMap = Ydata - (P(1) * Ddata + P(2));
w = BldMap > 0;
W = [w w ~w ~w]; A = [Ddata ones(nMeas, 1)]; A = [A A];
A_ = W .* A;
X = (A_' * A_) \ A_' * Ydata;
R.Z1 = [w ~w];
R.Alpha = [X(1), X(3)];
R.Beta = [X(2), X(4)];

% Generate collinear obstacles set
obstacles = (1:lenX*lenY)';
nObst = length(obstacles);
cols = cell(nObst, 2);   % collinear measurements of each obstacle
covB = cell(nMeas, 2);
for id = 1:nMeas
    DronePosMeter = Xdata(id, 1:3);
    DronePosPixel = [floor(DronePosMeter(1:2)/meterPerPixel)+1, DronePosMeter(3)];
    UserPosMeter = Xdata(id, 4:6);
    UserPosPixel = [floor(UserPosMeter(1:2)/meterPerPixel)+1, UserPosMeter(3)];
    [covBlds, covBldZs] = covBldZ(DronePosPixel, UserPosPixel, lenX);
    covBlds = covBlds(covBldZs > 0 & covBldZs <= droneHeightMap);
    covBldZs = covBldZs(covBldZs > 0 & covBldZs <= droneHeightMap);
    covB{id, 1} = covBlds;
    covB{id, 2} = covBldZs;
    for i = 1:length(covBlds)
        j = covBlds(i);
        ib = find(obstacles == j, 1, 'first');
        if ~isempty(ib)
            if ~isempty(cols{ib, 1})
                cols{ib, 1} = [cols{ib, 1} id];
                cols{ib, 2} = [cols{ib, 2} covBldZs(i)];
            else
                cols{ib, 1} = id;
                cols{ib, 2} = covBldZs(i);
            end
        end
    end
end

%% LOCAL POLY RECONSTRUCTION ----------------------------------------------
% I perform building height estimation on a multiresolution map, just to
% accelerate the algorithm. Building height estimation on coarse
% (discretized) map can be a good intialization for a fine map.
MAXLOOP = 9;
tol = 1e-3;
metrics = cell(MAXLOOP, 6);
Maps.droneHeightMap = droneHeightMap;
Maps.neighbourMat = neighbourMat;

% -- coarse resolution --
dsfactor = 4;       % Downsampling factor
[smallMaps, S, Hs] = downsampleMaps(Maps, dsfactor, R);
R.Hs = Hs; R.S = S;
sBldMapHat = ones(size(smallMaps.BldMapZ)) * min(droneHeightMap, ...
    max(droneHeightSample));

for i = 1:MAXLOOP
    [sBldMapHat, w] = virObsStep(R, smallMaps, sBldMapHat);
    W = [w w 1-w 1-w];
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    if norm([X(1) X(3) X(2) X(4)] - [R.Alpha R.Beta]) < tol
        break
    end
    R.Z1 = [-abs(Ydata-(X(1)*Ddata+X(2))) -abs(Ydata-(X(3)*Ddata+X(4)))];
    R.Alpha = [X(1), X(3)];
    R.Beta = [X(2), X(4)];
end

% -- medium resolution --
dsfactor2 = 2;  
[mediumMaps, S, Hs] = downsampleMaps(Maps, dsfactor2, R);
R.Hs = Hs; R.S = S;
mBldMapHat = upsamplematrix(sBldMapHat, size(mediumMaps.BldPosMat));
mBldMapHat(mediumMaps.BldPosMat < 0.1) = 0;

for i = 1:MAXLOOP
    [mBldMapHat, w] = virObsStep(R, mediumMaps, mBldMapHat);
    W = [w w 1-w 1-w];
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    if norm([X(1) X(3) X(2) X(4)] - [R.Alpha R.Beta]) < tol
        break
    end
    R.Z1 = [-abs(Ydata-(X(1)*Ddata+X(2))) -abs(Ydata-(X(3)*Ddata+X(4)))];
    R.Alpha = [X(1), X(3)];
    R.Beta = [X(2), X(4)];
end

% -- fine resolution (3 meter) --
dsfactor3 = 1;
[Maps, S, Hs] = downsampleMaps(Maps, dsfactor3, R);
R.Hs = Hs; R.S = S;
BldMapHat = upsamplematrix(mBldMapHat, size(Maps.BldPosMat));
BldMapHat(Maps.BldPosMat < 0.2) = 0;

h  = waitbar(0, 'Estimate building heights');
for i = 1:MAXLOOP
    waitbar(i / MAXLOOP, h);
    W = [w w 1-w 1-w];
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    % Parameters denoise
    D1 = abs(Ydata - (X(1)*Ddata+X(2))); w1 = (D1 < std(D1(w>0.5))*4)&(w>0.5);
    D2 = abs(Ydata - (X(3)*Ddata+X(4))); w2 = (D2 < std(D2(w<0.5))*4)&(w<0.5);
    W = [w1 w1 w2 w2]; A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    metrics{i, 1} = X; metrics{i, 2} = norm(A_*X-Ydata)^2;
    if norm([X(1) X(3) X(2) X(4)] - [R.Alpha R.Beta]) < tol && ...
            i ~= 1 && norm(metrics{i, 2} - metrics{i-1, 2}) < tol
        break
    end
    R.Z1 = [-abs(Ydata-(X(1)*Ddata+X(2))) -abs(Ydata-(X(3)*Ddata+X(4)))];
    R.Alpha = [X(1), X(3)];
    R.Beta = [X(2), X(4)];
    [BldMapHat, w] = virObsStep(R, Maps, BldMapHat);
end
close(h);

% hf = showmap(BldMapHat + FolMapHat, Maps.meterPerPixel, 8);
% title('Esitmated building and foliage height');

if kmeansOn == true
    c_i = D1 ./ abs((X(1)*Ddata+X(2)) - (X(3)*Ddata+X(4)));
    obs_clu = zeros(nObst, droneHeightMap);
    for i = 1:droneHeightMap
        for j = 1:nObst
            obs_clu(j, i) = mean(c_i(cols{j, 1}(ceil(cols{j, 2}) == i)));
        end
    end
    obs_clu(isnan(obs_clu)) = 0;
    BldInd = zeros(lenX, lenY);
    BldInd(Maps.BldPosMat > 0) = kmeans(obs_clu, 5);
end

figure
subplot(2,2,1); hold on
plot(Ddata(w1), Ydata(w1), 'k.')
plot(Ddata(w2), Ydata(w2), 'b.')
subplot(2,2,2); hold on
plot(Ddata(w2), Ydata(w2), 'k.')
plot(Ddata(w1), Ydata(w1), 'b.')
subplot(2,2,[3,4]); hold on
plot(Ddata, Ydata, '.')

heights = BldMapHat;
for i = 1:length(userPosAll)
    x = round(userPosAll(i, 1) / meterPerPixel);
    y = round(userPosAll(i, 2) / meterPerPixel);
    if x == 0 || y == 0, continue, end
    heights(x, y) = 0;
end
% heights(heights >= droneHeightMap | heights == droneHeightMap/2) = 0;
% indicator = heights >= UserHeight;
% heights = medfilt1(heights .* (heights <= DroneHeight)) .* indicator;
figure; showmap(heights, meterPerPixel);
title(sprintf('Height from L model at %d meters', round(droneHeightMap)));

DroneNum = DroneNum * length(droneHeightSample);

%% RADIOMAP AND PERFORMANCE -----------------------------------------------
% Channel model 
C.A1 = X(1); C.B1 = X(2); C.S1 = 0;
C.A2 = 0; C.B2 = 0; C.S2 = 0; 
C.A3 = X(3); C.B3 = X(4); C.S3 = 0;
ChannelModel.C = C;
ChannelModel.los_nlos_trans = 0;
ChannelModel.noise = 0;

metrics_mse = []; metrics_mae = [];
XYkr = [vec(meshgrid(1:lenX, 1:lenY)'), vec(meshgrid(1:lenY, 1:lenX)), ...
    ones(lenX * lenY, 1) * droneHeightMap / meterPerPixel];
Xann = [vec(meshgrid(1:lenX, 1:lenY)'), vec(meshgrid(1:lenY, 1:lenX)) ...
    * meterPerPixel, ones(lenX * lenY, 1) * droneHeightMap];
net = fitnet([32 32 16 16 8 8 4 4]); %fitnet([(4:8).^2 (8:-1:2).^2]);
net = train(net, [Xdata Ddata X(1)*Ddata+X(2) X(3)*Ddata+X(4)]', Ydata', ...
    'useParallel', 'no', 'showResources', 'no');
j = 1;
for i = 1:length(PosUE)
    if userPosAll(i, 1) == 0 || userPosAll(i, 2) == 0, continue, end
    if ~ismember(i, testPositions), j = j + 1; continue, end
    Gtru = Gtrue{i}';
    % Radiomap reconstruction
    Xdata_i = [Xdata((j-1)*DroneNum+1:j*DroneNum,1:2)/meterPerPixel, ...
        Xdata((j-1)*DroneNum+1:j*DroneNum,3)/meterPerPixel];
    Ydata_i = Ydata((j-1)*DroneNum+1:j*DroneNum);
    idx = randperm(length(Xdata_i), round(length(Xdata_i) / sampleRate));
    Xdata_i = Xdata_i(idx, :); Ydata_i = Ydata_i(idx, :);
    d_variogram = variogram(Xdata_i(:, 1:2), Ydata_i);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    if nMeas > 1e9
        idx = randperm(length(Xdata_i), round(length(Xdata_i) / sampleRate));
        Gkri = kriging(vstruct,Xdata_i(idx,:),false,Ydata_i(idx),XYkr,false);
    else
        idx = randperm(length(Xdata), round(length(Xdata) / sampleRate));
        Gkri = kriging(vstruct, Xdata(idx,:)/meterPerPixel, false, Ydata(idx), ...
        [XYkr ones(length(XYkr), 3).*userPosAll(i, :)/meterPerPixel], false);
    end
    Xann_upos = ones(length(XYkr), 3).*userPosAll(i, :);
    Xann_dist = log10(vecnorm(Xann - Xann_upos, 2, 2));
    Xann_features = [Xann Xann_upos ...
        Xann_dist X(1)*Xann_dist+X(2) X(3)*Xann_dist+X(4)]';
    Gann = net(Xann_features);
    Gann = reshape(Gann, [lenX lenY]);
    Gkri = reshape(Gkri, [lenX lenY]);
    Gknn = powerMapReconKNN(R,Ydata,userPosAll(i, :),droneHeightMap,Maps);
    Gbld = powerMapReconBld(R,heights,heights,userPosAll(i, :),droneHeightMap,Maps);
    Gbdn = Gbld;
    % Radiomap reconstruction with residual
    if residualOn == true
    [~, id] = min(abs(droneHeightSample - droneHeightMap));
    droneHeightRes = droneHeightSample(id)/meterPerPixel;
    ResX = floor(Xdata_i(:, 1:2)) + 1;
    ResX = ResX(Xdata_i(:, 3) == droneHeightRes, :);
    ResG = Ydata_i(Xdata_i(:, 3) == droneHeightRes);
    Gbld = powerMapReconBld(R,heights,heights,userPosAll(i, :),droneHeightMap,Maps);
    Gbld = imgaussfilt(Gbld, 1);
    for k = 1:length(ResX)
        ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
    end
    d_variogram = variogram(ResX, ResG);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    Rkri = kriging(vstruct, ResX(:,1), ResX(:,2), ...
        ResG, meshgrid(1:lenX, 1:lenY)', meshgrid(1:lenY, 1:lenX));
    Gbld = Gbld + reshape(Rkri, [lenX lenY]);
    end
    if residsecOn == true
    ResX = Xdata / meterPerPixel;
    P = (X(1)*Ddata + X(2) - Ydata)./(X(1)*Ddata + X(2) - X(3)*Ddata - X(4));
    P(P < 0) = 0; P(P > 1) = 1;
    ResG = Ydata - w1.*(X(1)*Ddata + X(2)) - w2.*(X(3)*Ddata + X(4));
    d_variogram = variogram(Xdata((j-1)*DroneNum+1:j*DroneNum,1:2)/meterPerPixel, ...
        Ydata((j-1)*DroneNum+1:j*DroneNum));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    ResRec = kriging(vstruct, ResX, false, ResG, ...
    [XYkr ones(length(XYkr), 3).*userPosAll(i, :)/meterPerPixel], false);
    ResRec = reshape(ResRec, [lenX lenY]);
    Gbld = Gbld + ResRec;
    end
    if residthiOn == true
    Rind = true(length(ResX), 1);
    for k = 1:length(ResX)
        ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
        Rind(k) = Gind(ResX(k, 1), ResX(k, 2));
    end
    d_variogram = variogram(ResX(Rind,:), ResG(Rind));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    LOSSeg = kriging(vstruct, ResX(Rind,1), ResX(Rind,2), ...
        ResG(Rind), meshgrid(1:lenX, 1:lenY)', meshgrid(1:lenY, 1:lenX));
    d_variogram = variogram(ResX(~Rind,:), ResG(~Rind));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    NLOSSeg = kriging(vstruct, ResX(~Rind,1), ResX(~Rind,2), ...
        ResG(~Rind), meshgrid(1:lenX, 1:lenY)', meshgrid(1:lenY, 1:lenX));
    Gbld = Gbld + LOSSeg .* Gind + NLOSSeg .* ~Gind;
    end
    if residfouOn == true
    ResY = []; ResR = [];
    for id = 1:length(droneHeightSample)
        h = droneHeightSample(id)/meterPerPixel;
        ResX = [floor(Xdata_i(:, 1:2)) + 1 Xdata_i(:, 3)];
        ResX = ResX(Xdata_i(:, 3) == h, :);
        ResG = Ydata_i(Xdata_i(:, 3) == h);
        Gbld = powerMapReconBld(R,heights,heights,userPosAll(i, :),h,Maps);
        Gbld = imgaussfilt(Gbld, 3);
        for k = 1:length(ResX)
            ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
        end
        x = ResX(:, 1) - userPosAll(i, 1)/meterPerPixel;
        y = ResX(:, 2) - userPosAll(i, 2)/meterPerPixel;
        h = h/meterPerPixel - userPosAll(i, 3)/meterPerPixel;
        l = vecnorm([x y], 2, 2);
        theta1 = atan2(y,x); theta2 = atan2(h,l);
        ResR = [ResR; min(theta1+pi,pi-theta1)/pi*180 ...
            min(theta2+pi,pi-theta2)/pi*180 sqrt(l.^2 + h^2)];
        ResY = [ResY; ResG];
    end
    d_variogram = variogram(ResR, ResY);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    Rkri = krigingR(vstruct, Xdata((j-1)*DroneNum+1:j*DroneNum,1:3)/ ...
        meterPerPixel, ResY, XYkr, userPosAll(i, :)/meterPerPixel);
    Gbld = powerMapReconBld(R,heights,heights,userPosAll(i, :),droneHeightMap,Maps);
    Gbld = imgaussfilt(Gbld, 3) + reshape(Rkri, [lenX lenY]);
    end
    if residfivOn == true
    [~, id] = min(abs(droneHeightSample - droneHeightMap));
    droneHeightRes = droneHeightSample(id)/meterPerPixel;
    ResX = floor(Xdata_i(:, 1:2)) + 1;
    ResX = ResX(Xdata_i(:, 3) == droneHeightRes, :);
    ResG = Ydata_i(Xdata_i(:, 3) == droneHeightRes);
    Gind = powMapRecBldSof(R,heights,[],userPosAll(i, :),droneHeightMap,Maps,[]);
    Gind = imgaussfilt(Gind, 3);
    Glos = powMapRecBldSof(R,heights-realmax,[],userPosAll(i, :),droneHeightMap,Maps);
    Gnlo = powMapRecBldSof(R,heights+realmax,[],userPosAll(i, :),droneHeightMap,Maps);
    Gbld = Gind .* Glos + (1 - Gind) .* Gnlo;
    Gbdn = Gbld;
    for k = 1:length(ResX)
        ResG(k) = ResG(k) - Gbld(ResX(k, 1), ResX(k, 2));
    end
    d_variogram = variogram(ResX, ResG);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    Rkri = kriging(vstruct, ResX(:,1), ResX(:,2), ...
        ResG, meshgrid(1:lenX, 1:lenY)', meshgrid(1:lenY, 1:lenX));
    Gbld = Gbld + reshape(Rkri, [lenX lenY]);
    end
    % Radiomap reconstruction with segmentation
    if segmentOn == true
    d_variogram = variogram(Xdata(w1, :)/meterPerPixel, Ydata(w1));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    LOSSeg = kriging(vstruct,Xdata(w1, :)/meterPerPixel, false, Ydata(w1), ...
        [XYkr ones(length(XYkr), 3).*userPosAll(i, :)/meterPerPixel], false);
    d_variogram = variogram(Xdata(w2, :)/meterPerPixel, Ydata(w2));
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    NLOSSeg = kriging(vstruct,Xdata(w2, :)/meterPerPixel, false, Ydata(w2), ...
        [XYkr ones(length(XYkr), 3).*userPosAll(i, :)/meterPerPixel], false);
    LOSSeg = reshape(LOSSeg, [lenX lenY]);
    NLOSSeg = reshape(NLOSSeg, [lenX lenY]);
    Gind = powerMapReconBld(R,heights,heights,userPosAll(i, :),droneHeightMap,Maps,1);
    Gbld = LOSSeg .* Gind + NLOSSeg .* ~Gind;
    end
    % Radiomap reconstruction with residual and segmentation
    if segresidOn == true
    Xlos = Xdata(w1, :) / meterPerPixel; 
    Ylos = Ydata(w1) - Ddata(w1) * X(1) - X(2);
    Xnlos = Xdata(w2, :) / meterPerPixel;
    Ynlos = Ydata(w2) - Ddata(w2) * X(3) - X(4);
    d_variogram = variogram(Xlos, Ylos);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    LOSSeg = kriging(vstruct, Xlos, false, Ylos, ...
        [XYkr ones(length(XYkr), 3).*userPosAll(i, :)/meterPerPixel], false);
    d_variogram = variogram(Xnlos, Ynlos);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
        [], [], [], 'model', 'exponential', 'plotit', false);
    NLOSSeg = kriging(vstruct, Xnlos, false, Ynlos, ...
        [XYkr ones(length(XYkr), 3).*userPosAll(i, :)/meterPerPixel], false);
    LOSSeg = reshape(LOSSeg, [lenX lenY]);
    NLOSSeg = reshape(NLOSSeg, [lenX lenY]);
    Gind = powerMapReconBld(R,heights,heights,userPosAll(i, :),droneHeightMap,Maps,1);
    Gbld = Gbld + LOSSeg .* Gind + NLOSSeg .* ~Gind;
    end
    % Performance and figures
    xm = length(Zmat);
    if xm <= lenX
        Gbld = Gbld(round(1:lenX/xm:lenX),round(1:lenX/xm:lenX));
        Gkri = Gkri(round(1:lenX/xm:lenX),round(1:lenX/xm:lenX));
        Gknn = Gknn(round(1:lenX/xm:lenX),round(1:lenX/xm:lenX));
        Gann = Gann(round(1:lenX/xm:lenX),round(1:lenX/xm:lenX));
        Gbdn = Gbdn(round(1:lenX/xm:lenX),round(1:lenX/xm:lenX));
    else
        Gtru = Gtru(round(1:xm/lenX:xm),round(1:xm/lenX:xm));
    end
    if plotOn == true
    figure(300 + i);
    subplot(2, 3, 1); showmap(Gkri, meterPerPixel);
    subplot(2, 3, 2); showmap(Gknn, meterPerPixel);
    subplot(2, 3, 3); showmap(Gann, meterPerPixel);
    subplot(2, 3, 4); showmap(Gbdn, meterPerPixel);
    subplot(2, 3, 5); showmap(Gbld, meterPerPixel);
    subplot(2, 3, 6); showmap(Gtru, meterPerPixel);
    % figure(400 + i); showmap(reshape(Rkri, [lenX lenY]), 9);
    input('Position '+string(i)+': '+string(mean(abs(vec(Gtru-Gkri)))) ...
        +' vs. '+string(mean(abs(vec(Gtru-Gknn)))) ...
        +' vs. '+string(mean(abs(vec(Gtru-Gbld)))) ...
        +'; '+string(mse(Gtru-Gkri)) +' vs. '+string(mse(Gtru-Gknn)) ...
        +' vs. '+string(mse(Gtru-Gbld)));
    end
    metrics_mse = [metrics_mse; mse(Gtru-Gkri) mse(Gtru-Gknn) ...
        mse(Gtru-Gann) mse(Gtru-Gbdn) mse(Gtru-Gbld)];
    metrics_mae = [metrics_mae; mean(abs(vec(Gtru-Gkri))) ...
        mean(abs(vec(Gtru-Gknn))) mean(abs(vec(Gtru-Gann))) ...
        mean(abs(vec(Gtru-Gbdn))) mean(abs(vec(Gtru-Gbld)))];
    j = j + 1;
end
disp([mean(metrics_mae, 1) mean(metrics_mse, 1) nUser]);
end
%% LOCALIZATION APPLICATION -----------------------------------------------
input('LOCALIZATION?');
for nDroneLoc = 6:20
plotOn = false;
eps_h = 1e-6;
learningRate = 1e-1;
MAXLOOP = 5e4;

Maps.BldMapZ = heights;
Maps.BldPosMat = ones(lenX, lenY);
Maps.FolPosMat = zeros(lenX, lenY);
Maps.meterPerPixel = meterPerPixel;
C.A1 = X(1); C.B1 = X(2); C.S1 = 0;
C.A2 = 0; C.B2 = 0; C.S2 = 0; 
C.A3 = X(3); C.B3 = X(4); C.S3 = 0;
ChannelModel.C = C;
ChannelModel.los_nlos_trans = 0;
ChannelModel.noise = 0;
DronePosArray = [floor(Xdata(:, 1:2)/meterPerPixel)+1 Xdata(:, 3)];
UserPosArray = [floor(Xdata(:, 4:5)/meterPerPixel)+1 Xdata(:, 6)];
Gainrow = powersample([DronePosArray UserPosArray], Maps, ChannelModel);
Rdata = Gainrow' - Ydata;
d_variogram = variogram(Xdata, Rdata);
[~, ~, ~, vstruct] = variogramfit(d_variogram.distance, d_variogram.val, ...
    [], [], [], 'model', 'exponential', 'plotit', false);

testUsers = setdiff(1:length(PosUE), userids);
testPhaseSize = 3;
testUserSize = length(testUsers);
metrics_err = zeros(testPhaseSize, testUserSize);
metrics_obj = zeros(testPhaseSize, testUserSize);
userPosAllXY = userPosAll(userids, 1:2);
for testPhaseNum = 1:testPhaseSize
parfor u = 1:testUserSize
    userId = testUsers(u);
    droneHeight = 50;
    x = PosUE(userId,1); y = PosUE(userId,2); z = PosUE(userId,3);
    I = (DATA.RadioMap(:, 6) == droneHeight) ...
        & (DATA.RadioMap(:, 1) == x) ...
        & (DATA.RadioMap(:, 2) == y)...
        & (DATA.RadioMap(:, 3) == z);
    Rm2D = DATA.RadioMap(I, :);
    Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
    [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
    Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
    UserPos = [x-min(Xvec)+1 y-min(Yvec)+1 z];
    Xmat = Xmat - min(Xvec) + 1;
    Ymat = Ymat - min(Yvec) + 1;
    % Data pre-processing
    Zmat(isnan(Zmat)) = min(min(Zmat));
    Xvec = Xmat(:);
    Yvec = Ymat(:);
    Zvec = Zmat(:);
    %idx = 1:ceil(length(Xvec)/nDroneLoc):length(Xvec);
    idx = sort(randperm(length(Xvec), nDroneLoc));
    Xvec = Xvec(idx);
    Yvec = Yvec(idx);
    Zvec = Zvec(idx);
    Hvec = ones(length(idx), 1) * droneHeight;
    % Initialize estimation of user's position
    losVar = 4 * std(D1(w1));
    nlosVar = 4 * std(D2(w2));
    d = zeros(length(idx), 1);
    w3 = true(length(idx), 1);
    w4 = true(length(idx), 1);
    metrics_xupos = [];
    [~, id] = max(Zvec);
    xuPos = [Xvec(id) Yvec(id)]+rand*30;%UserPos(1:2)+rand*30;%
    cnt = 0;
    xuPos0 = xuPos * 2;
    while cnt < MAXLOOP && norm(xuPos0 - xuPos, 'inf') >= eps_h
        xuPos0 = xuPos;
        objGrad = 0;
        Zvh = Zvec;
        % First order derivative
        g_los = @(x0,x,i) -(2*R.Alpha(1)*(R.Beta(1)+(R.Alpha(1)*log(sqrt((z- ...
            droneHeight)^2+(x-x0)*(x-x0)')))/log(10)-Zvh(i)))/(((z- ...
            droneHeight)^2+(x-x0)*(x-x0)')*log(10))*(x-x0);
        g_nlos = @(x0,x,i) -(2*R.Alpha(2)*(R.Beta(2)+(R.Alpha(2)*log(sqrt((z- ...
            droneHeight)^2+(x-x0)*(x-x0)')))/log(10)-Zvh(i)))/(((z- ...
            droneHeight)^2+(x-x0)*(x-x0)')*log(10))*(x-x0);
        for i = 1:length(idx)
            DronePos = [floor([Xvec(i) Yvec(i)]/meterPerPixel)+1, droneHeight];
            UE1Pos = [floor(xuPos/meterPerPixel)+1, z];
            noBld = true;
            [covBlds, covBldZs] = covBldZ(DronePos, UE1Pos, lenX);
            covBlds = covBlds(covBldZs > z & covBldZs <= droneHeightMap);
            covBldZs = covBldZs(covBldZs > z & covBldZs <= droneHeightMap);
            for j = 1:length(covBlds)
                k = covBlds(j);
                yb = max(min(floor((k - 1) / lenX) + 1, lenY), 1);
                xb = max(min(k - (yb - 1) * lenX, lenX), 1);
                zj = covBldZs(j);
                if zj < heights(xb, yb)  
                    noBld = false;
                end
            end
            d(i) = log10(norm([Xvec(i) Yvec(i) droneHeight] - [xuPos z], 2));
            noBld = noBld & (((X(1)*d(i)+X(2)) - Zvh(i)) < losVar);
            w4(i) = ~noBld & (Zvh(i) - ((X(3)*d(i)+X(4))) < nlosVar);
            w3(i) = ~w4(i) & (((X(1)*d(i)+X(2)) - Zvh(i)) < losVar);
            if w3(i)
                objGrad = objGrad + g_los(xuPos, [Xvec(i) Yvec(i)], i);
            elseif w4(i)
                objGrad = objGrad + g_nlos(xuPos, [Xvec(i) Yvec(i)], i);
            end
        end
        % Normalized gradient
        objGrad = objGrad / sum(w3 | w4);
        % Forced objective value descent
        f_obj = @(x0) norm(R.Beta(1)+(R.Alpha(1)*log10(sqrt((z-droneHeight)^2 ...
            +vecnorm([Xvec(w3) Yvec(w3)]-x0,2,2).^2)))-Zvec(w3)) + ...
            norm(R.Beta(2)+(R.Alpha(2)*log10(sqrt((z-droneHeight)^2 ...
            +vecnorm([Xvec(w4) Yvec(w4)]-x0,2,2).^2)))-Zvec(w4));
        xuPos = xuPos - learningRate * objGrad;
        xuPos(xuPos > lenX * meterPerPixel) = lenX * meterPerPixel;
        % Metrics and figures
        metrics_xupos = [metrics_xupos; xuPos];
        cnt = cnt + 1;
        if cnt > 1e3 && (norm(metrics_xupos(cnt - 99, :) - xuPos, 'inf') <= 1e-3 ...
                || norm(metrics_xupos(cnt - 98, :) - xuPos, 'inf') <= 1e-3 ...
                || norm(mean(metrics_xupos(cnt - 999:end - 500, :), 1) - ...
                mean(metrics_xupos(cnt - 500:end, :), 1), 'inf') <= 1e-3)
            break
        end
        if mod(cnt, 1e4) == 1 && plotOn
            figure
            subplot(2,2,1); hold on
            plot(d(w3), Zvec(w3), 'k.')
            plot(d(w4), Zvec(w4), 'r.')
            subplot(2,2,2); hold on
            plot(d(w4), Zvec(w4), 'k.')
            plot(d(w3), Zvec(w3), 'r.')
            subplot(2,2,[3,4]); hold on
            plot(d, Zvec, '.')
            pause(0.5);
        end
    end

    lowBound = max(xuPos - 50, ones(1, 2));
    upBound = min(xuPos + 50, size(heights) * meterPerPixel);
    minObjVal = f_obj(xuPos);
    for xuPosX = lowBound(1):10:upBound(1)
        for xuPosY = lowBound(2):10:upBound(2)
            xuPos_bak = xuPos;
            xuPos = [xuPosX xuPosY];
            [~, id] = min(vecnorm(userPosAllXY - xuPos, 2, 2));
            ids = (sum(Xdata(:, 4:5) == userPosAllXY(id, :), 2) == 2);
            XdataL = Xdata(ids, :); RdataL = Rdata(ids, :);
            Zvh = Zvec + kriging(vstruct, XdataL, false, RdataL, ...
                [Xvec Yvec Hvec [xuPos z] .* ones(length(idx), 3)], false);
            for i = 1:length(idx)
                DronePos = [floor([Xvec(i) Yvec(i)]/meterPerPixel)+1, droneHeight];
                UE1Pos = [floor(xuPos/meterPerPixel)+1, z];
                noBld = true;
                [covBlds, covBldZs] = covBldZ(DronePos, UE1Pos, lenX);
                covBlds = covBlds(covBldZs > z & covBldZs <= droneHeightMap);
                covBldZs = covBldZs(covBldZs > z & covBldZs <= droneHeightMap);
                for j = 1:length(covBlds)
                    k = covBlds(j);
                    yb = max(min(floor((k - 1) / lenX) + 1, lenY), 1);
                    xb = max(min(k - (yb - 1) * lenX, lenX), 1);
                    zj = covBldZs(j);
                    if zj < heights(xb, yb)  
                        noBld = false;
                    end
                end
                d(i) = log10(norm([Xvec(i) Yvec(i) droneHeight] - [xuPos z], 2));
                noBld = noBld & (((X(1)*d(i)+X(2)) - Zvh(i)) < losVar);
                w4(i) = ~noBld & (Zvh(i) - ((X(3)*d(i)+X(4))) < nlosVar);
                w3(i) = ~w4(i) & (((X(1)*d(i)+X(2)) - Zvh(i)) < losVar);
            end
            f_obj = @(x0) ( norm(R.Beta(1)+(R.Alpha(1)*log10(sqrt((z- ...
                droneHeight)^2+vecnorm([Xvec(w3) Yvec(w3)]-x0,2,2).^2))) ...
                -Zvh(w3)) + norm(R.Beta(2)+(R.Alpha(2)*log10(sqrt((z- ...
                droneHeight)^2+vecnorm([Xvec(w4) Yvec(w4)]-x0,2,2).^2))) ...
                -Zvh(w4)) ) / sum(w3 | w4);
            objVal = f_obj(xuPos);
            xuPos = xuPos_bak;
            if minObjVal > objVal
                xuPos = [xuPosX xuPosY];
                minObjVal = objVal;
            end
        end
    end
    metrics_err(testPhaseNum, u) = norm(xuPos - UserPos(1:2));
    metrics_obj(testPhaseNum, u) = minObjVal;
end
end
[~, id] = min(metrics_obj, [], 1);
metrics_err = metrics_err(id + (0:testPhaseNum:testPhaseNum*(testUserSize-1)));
disp('MAE for '+string(nDroneLoc)+' samples: '+string(mean(metrics_err)))
end
%% BASELINE OF APPLICATION ------------------------------------------------
for nDroneLoc = 6:20%10:10:100
eps_h = 1e-6;
learningRate = 1e-1;
MAXLOOP = 5e4;
P = polyfit(Ddata, Ydata, 1);
metrics_xerr_bas = [];
for userId = 1:100
    droneHeight = 50;
    x = PosUE(userId,1); y = PosUE(userId,2); z = PosUE(userId,3);
    I = (DATA.RadioMap(:, 6) == droneHeight) ...
        & (DATA.RadioMap(:, 1) == x) ...
        & (DATA.RadioMap(:, 2) == y)...
        & (DATA.RadioMap(:, 3) == z);
    Rm2D = DATA.RadioMap(I, :);
    Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
    [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
    Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
    UserPos = [x-min(Xvec)+1 y-min(Yvec)+1 z];
    Xmat = Xmat - min(Xvec) + 1;
    Ymat = Ymat - min(Yvec) + 1;
    % Data pre-processing
    Zmat(isnan(Zmat)) = min(min(Zmat));
    Xvec = Xmat(:);
    Yvec = Ymat(:);
    Zvec = Zmat(:);
    %idx = 1:ceil(length(Xvec)/nDroneLoc):length(Xvec);
    idx = sort(randperm(length(Xvec), nDroneLoc));
    Xvec = Xvec(idx);
    Yvec = Yvec(idx);
    Zvec = Zvec(idx);
    
    [~, id] = max(Zvec);
    xuPos = [Xvec(id) Yvec(id)];
    sv_pos = [Xvec Yvec droneHeight * ones(length(idx), 1)]';
    pr = 10 .^ ((Zvec - P(2)) / P(1)); pr = pr';
%{
    cnt = 0;
    xuPos0 = xuPos * 2;
    while cnt < MAXLOOP && norm(xuPos0 - xuPos, 'inf') >= eps_h
        xuPos0 = xuPos;
        g_obj = @(x0,x,di) -(2*(((50-z)^2+(x-x0)*(x-x0)').^0.5-di))/ ...
            ((50-z)^2+(x-x0)*(x-x0)').^0.5*(x-x0);
        objGrad = 0
        for i = 1:length(idx)
            objGrad = objGrad + g_obj(xuPos, [Xvec(i) Yvec(i)], pr(i));
        end
        objGrad = objGrad / length(idx);
        xuPos = xuPos - learningRate * objGrad
        cnt = cnt + 1;
    end
%}
    [N1, ~] = Trilateration(sv_pos, pr, diag(ones(1,length(pr))));
    xuPos = N1(2:3, 1)';
    xuPos(xuPos > lenX * meterPerPixel) = lenX * meterPerPixel;
    xuPos(xuPos < 0) = 0;
    %disp([xuPos Xvec(id) Yvec(id) UserPos(1:2) norm(xuPos - UserPos(1:2))]);
    metrics_xerr_bas = [metrics_xerr_bas; norm(xuPos - UserPos(1:2))];
end
disp('MAE for '+string(nDroneLoc)+' samples: '+string(mean(metrics_xerr_bas)))
end
