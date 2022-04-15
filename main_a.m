% Air-to-ground power map reconstruction via virtual building height
% estimation and terrain classificaiton. 
%
% MODIFICATION HISTORY:
% Feature 1: 3 propagation classification, LOS, obstructed-LOS, and NLOS
% (corresponding to free space, foliage blockage, and building blockage)
%
% Feature 1.1: Uniform random noise between propagation segments
%
% Feature 2: Multiple resolution techniques to accelarate the learning
% algorithm
%
% Work in progress for the WI-UAV paper workshop paper. 
% Author: Junting Chen (juntingc@usc.edu)
% Date: June 28, 2017. 

clear
close all
addpath Sub
metrics_mae = zeros(1, 4);
for testPhaseNum = 1:size(metrics_mae, 1)

% Topology parameters
Nue = 50;       % Default 100
Ndrone = 200;    % Default 400
UE1 = [190, 123]; % meter
% Number of user/UAV for training
% 22.3615   26.5766   31.5864   37.5405   44.6170
% 53.0274   63.0231   74.9031   89.0225  100.0000

DroneHeightMap = 50;        % Drone height for air-to-ground power map (Reconstruction stage)
DroneHeightSample = 50;     % Drone height for learning (Learning stage)

noise = 9;
% Channel model 
C.A1 = - 22; % LOS (LTE TR36.814)
C.B1 = - 28;
C.S1 = noise;    % Model variance / shadowing in RMS (std)

C.A2 = - 28; % Obstructed LOS (foliage etc.)
C.B2 = - 24;
C.S2 = 3 * noise; 

C.A3 = - 36; % NLOS (LTE TR36.814)
C.B3 = - 22;
C.S3 = noise;

los_nlos_trans = 3;         % Size of transition region between segments
                            % This is to model the "edge effect", the
                            % signal should have a "smooth" transition from
                            % LOS to NLOS
ChannelModel.C = C;
ChannelModel.los_nlos_trans = los_nlos_trans;

% Urban topology data
DATA = load('BldMap_trueHeight.mat');
smallBldMap = ordfilt2(DATA.smallBldMap, 1, ones(3, 3));
meterPerPixel = DATA.dsMeterPerPixel * 5;
startp = 150;
stopp = 450;
BldMapZ = smallBldMap(startp:meterPerPixel:stopp, startp:meterPerPixel:stopp).';
BldPosMat = (BldMapZ > 6);
% DATA = load('BuildingMap_3m.mat');
% meterPerPixel = DATA.dsMeterPerPixel;
% BldMapZ = DATA.smallBldMap(1:end, 1:end).';
% BldPosMat = DATA.BldAreaIndicator(1:end, 1:end).';
% BldMapZ = BldMapZ .* BldPosMat;

% % Topology data for test
% % Disable filter in advance
% BldMapZ = zeros(94, 94);
% BldMapZ(70, 40:60) = 30;
% BldMapZ(50, 40:60) = 35;
% BldMapZ(25, :) = 10;
% BldPosMat = zeros(94, 94);
% BldPosMat(70, 40:60) = 1;
% BldPosMat(50, 40:60) = 1;

[lenX, lenY] = size(BldMapZ);   % lenX = # of rows, corresponding to # of pixels in x-axis
                                % lenY = # of columns, corresponding ot #
                                % of pixels in y-axis
                                
% Setting - buliding, foliage, user positions (Extract 2D street map)
BldPos = find(BldPosMat > 0);
nBldPos = length(BldPos);

foliageThreshold = 5;  % meter (artificially identify foliage area)
FoliagePos = find((BldMapZ > foliageThreshold) & (BldPosMat == 0));
nFoliagePos = length(FoliagePos);
FolPosMat = zeros(size(BldMapZ));
FolPosMat(FoliagePos) = 1;

UserPosAll = find((BldMapZ <= foliageThreshold) & (BldPosMat == 0));
nUserPos = length(UserPosAll);
UserPosMat = zeros(size(BldMapZ));
UserPosMat(UserPosAll) = 1;

hf = showmap(BldMapZ, meterPerPixel, 1);
title('Building height map');
% hf = showmap(BldPosMat, meterPerPixel, 2);
% title('Building location map');
% hf = showmap(FolPosMat, meterPerPixel, 3);
% title('Foliage location map');
% hf = showmap(UserPosMat, meterPerPixel, 4);
% title('User position map');

Maps.BldMapZ = BldMapZ;
Maps.BldPosMat = BldPosMat;
Maps.FolPosMat = FolPosMat;
Maps.meterPerPixel = meterPerPixel;

UE1_pixel = [1, 1] + floor(UE1 / meterPerPixel);
UE1Pos = [UE1_pixel, BldMapZ(UE1_pixel(1), UE1_pixel(2))];  % [pixel, pixel, meter]
UE1PosMeter = [UE1, UE1Pos(3)];

GroupNum = 9;
GroupNumSqrt = floor(sqrt(GroupNum));
EdgeGroupMap = zeros(lenX, lenY);
EdgeGroupMap(1:floor(lenX / GroupNumSqrt) * GroupNumSqrt, ...
             1:floor(lenY / GroupNumSqrt) * GroupNumSqrt) =...
    repelem(reshape(1:GroupNum, [GroupNumSqrt GroupNumSqrt]),...
    floor(lenX / GroupNumSqrt), floor(lenY / GroupNumSqrt));

%% TRUE POWER MAP FOR USER 1 ----------------------------------------------
ChannelModel.noise = 0;
[DX, DY] = meshgrid((1:lenX), (1:lenY));
DX = DX.';  DY = DY.';
DronePosArray = [DX(:) DY(:) DroneHeightMap * ones(length(DX(:)), 1)];
UserPosArray = repmat(UE1Pos, length(DX(:)), 1);
GainVec = powersample([DronePosArray UserPosArray], Maps, ChannelModel);
Gain = reshape(GainVec, lenX, lenY);

hf = showmap(Gain, meterPerPixel, 5);    
title(sprintf('True power map at drone height %d meters', round(DroneHeightMap)));
pause(0.2);

%% DATA COLLECTION --------------------------------------------------------
ChannelModel.noise = 1;
Nsam = Nue * Ndrone;
Xdata = zeros(Nsam, 6); 
Ydata = zeros(Nsam, 1);
Ddata = zeros(Nsam, 1);
cnt = 0;
droneX = round((1/(2 * ceil(sqrt(Ndrone))):1/ceil(sqrt(Ndrone)):1) * lenX);
droneY = round((1/(2 * ceil(sqrt(Ndrone))):1/ceil(sqrt(Ndrone)):1) * lenY);
[droneXmat, droneYmat] = meshgrid(droneX, droneY);
for idrone = 1:Ndrone
    DronePos = [droneXmat(idrone), droneYmat(idrone), DroneHeightSample];
    UserPoseIndex0 = randperm(nUserPos);
    UserPoseIndex = UserPosAll(UserPoseIndex0(1:Nue));
    for iuser = 1:Nue
        userIndex = UserPoseIndex(iuser);
        yu = floor((userIndex - 1) / lenX) + 1;
        xu = userIndex - (yu - 1) * lenX;
        zu = BldMapZ(xu, yu);
        UserPos = [xu, yu, zu];
        
        if EdgeGroupMap(DronePos(1:2)) ~= EdgeGroupMap(ceil(UserPos(1:2)))
            continue;
        end

        gain = powersample([DronePos, UserPos], Maps, ChannelModel);

        cnt = cnt + 1;
        DronePosMeter = [DronePos(1:2) * meterPerPixel, DronePos(3)];
        UserPosMeter = [UserPos(1:2) * meterPerPixel, UserPos(3)];
        Xdata(cnt, :) = [DronePosMeter, UserPosMeter];
        Ydata(cnt, :) = gain;
        Ddata(cnt, :) = log10(norm([DronePos(1:2) * meterPerPixel, ... 
            DronePos(3)] - [UserPos(1:2) * meterPerPixel, UserPos(3)]));
    end
end

Xdata = Xdata(1:cnt, :);
Ydata = Ydata(1:cnt, :);
Ddata = Ddata(1:cnt, :);

Nmeas = cnt;
R.X = Xdata;
P = polyfit(Ddata, Ydata, 1);
Z = Ydata - (P(1) * Ddata + P(2));
w = Z > 0;
W = [w w ~w ~w]; A = [Ddata ones(Nmeas, 1)]; A = [A A];
A_ = W .* A;
X = (A_' * A_) \ A_' * Ydata;
R.Z1 = [w ~w];
R.Alpha = [X(1), X(3)];
R.Beta = [X(2), X(4)];

%% LOCAL GRADIENT RECONSTRUCTION ------------------------------------------
tic
% I perform building height estimation on a multiresolution map, just to
% accelerate the algorithm. Building height estimation on coarse
% (discretized) map can be a good intialization for a fine map.
MAXLOOP = 60;
tol = 1e-3;
% -- coarse resolution --
dsfactor = 4;       % Downsampling factor
smallMaps = downsampleMaps(Maps, dsfactor);
sBldMapHat = ones(size(smallMaps.BldMapZ)) * DroneHeightMap;
sFolMapHat = ones(size(smallMaps.BldMapZ)) * DroneHeightMap;

for i = 1:MAXLOOP
    [sBldMapHat, sFolMapHat, w] = virtualBldMap(R, smallMaps, sBldMapHat + sFolMapHat);
    W = [w w ~w ~w];
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    if norm([X(1) X(3) X(2) X(4)] - [R.Alpha R.Beta]) < tol
        break
    end
%     D1 = exp((Ydata - (X(1) * Ddata + X(2))));
%     D2 = exp(((X(3) * Ddata + X(4)) - Ydata));
    R.Z1 = [-abs(Ydata-(X(1)*Ddata+X(2))) -abs(Ydata-(X(3)*Ddata+X(4)))];
    R.Alpha = [X(1), X(3)];
    R.Beta = [X(2), X(4)];
%     disp(X);
end

% -- medium resolution --
dsfactor2 = 2;  
mediumMaps = downsampleMaps(Maps, dsfactor2);
mBldMapHat = upsamplematrix(sBldMapHat, size(mediumMaps.BldPosMat));
mBldMapHat(mediumMaps.BldPosMat < 0.1) = 0;
mFolMapHat = upsamplematrix(sFolMapHat, size(mediumMaps.FolPosMat));
mFolMapHat(mediumMaps.FolPosMat < 0.1) = 0;

for i = 1:MAXLOOP
    [mBldMapHat, mFolMapHat, w] = virtualBldMap(R, mediumMaps, mBldMapHat + mFolMapHat);
    W = [w w ~w ~w];
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    if norm([X(1) X(3) X(2) X(4)] - [R.Alpha R.Beta]) < tol
        break
    end
%     D1 = exp((Ydata - (X(1) * Ddata + X(2))));
%     D2 = exp(((X(3) * Ddata + X(4)) - Ydata));
    R.Z1 = [-abs(Ydata-(X(1)*Ddata+X(2))) -abs(Ydata-(X(3)*Ddata+X(4)))];
    R.Alpha = [X(1), X(3)];
    R.Beta = [X(2), X(4)];
%     disp(X);
end

% -- fine resolution (3 meter) --
BldMapHat = upsamplematrix(mBldMapHat, size(Maps.BldPosMat));
BldMapHat(Maps.BldPosMat < 0.2) = 0;
FolMapHat = upsamplematrix(mFolMapHat, size(Maps.FolPosMat));
FolMapHat(Maps.FolPosMat < 0.2) = 0;

for i = 1:MAXLOOP
    [BldMapHat, FolMapHat, w] = virtualBldMap(R, Maps, BldMapHat + FolMapHat);
    W = [w w ~w ~w];
    A_ = W .* A;
    X = (A_' * A_) \ A_' * Ydata;
    if norm([X(1) X(3) X(2) X(4)] - [R.Alpha R.Beta]) < tol
        break
    end
%     D1 = exp((Ydata - (X(1) * Ddata + X(2))));
%     D2 = exp(((X(3) * Ddata + X(4)) - Ydata));
    R.Z1 = [-abs(Ydata-(X(1)*Ddata+X(2))) -abs(Ydata-(X(3)*Ddata+X(4)))];
    R.Alpha = [X(1), X(3)];
    R.Beta = [X(2), X(4)];
%     disp(X);
end
hf = showmap(BldMapHat + FolMapHat, Maps.meterPerPixel, 8);
title('Esitmated building and foliage height');
% save_figure(hf, fileHead, '15', {'fig','jpg'});
% hf = showmap(FolMapHat, Maps.meterPerPixel, 16); title('Estimated foliage height');
% save_figure(hf, fileHead, '16', {'fig','jpg'});

% Parameters denoise
D1 = abs(Ydata-(X(1)*Ddata+X(2)));
w1 = (D1 < std(D1)*2) .* w;
D2 = abs(Ydata-(X(3)*Ddata+X(4)));
w2 = (D2 < std(D2)*2) .* ~w;
W = [w1 w1 w2 w2];
A_ = W .* A;
X = (A_' * A_) \ A_' * Ydata;
R.Alpha = [X(1), X(3)];
R.Beta = [X(2), X(4)];

Gloc = powerMapReconBld(R, BldMapHat, FolMapHat, UE1PosMeter, DroneHeightMap, Maps);
hf = showmap(Gloc, Maps.meterPerPixel, 9);
title(sprintf('L reconstruction at drone height %d meters', round(DroneHeightMap)));
pause(0.2);

toc

figure
hold on
plot(Ddata(w >= 0.5), Ydata(w >= 0.5), 'b*')
plot(Ddata(w <= 0.5), Ydata(w <= 0.5), 'g*')
hold off
figure
plot(Ddata, Ydata, '*')
figure
hold on
plot(Ddata(w <= 0.5), Ydata(w <= 0.5), 'b*')
plot(Ddata(w >= 0.5), Ydata(w >= 0.5), 'g*')
hold off
heights_lg = BldMapHat + FolMapHat;

%% RADIOMAP AND PERFORMANCE -----------------------------------------------
input('Performance metrics?')
metrics_mse = []; metrics_mae = [];
userPosInd = UserPosAll(UserPoseIndex0(1:Nue));
userPosAll = zeros(Nue, 3);
for iuser = 1:Nue
    userIndex = userPosInd(iuser);
    yu = floor((userIndex - 1) / lenX) + 1;
    xu = userIndex - (yu - 1) * lenX;
    zu = BldMapZ(xu, yu);
    userPosAll(iuser, :) = [xu * meterPerPixel, yu * meterPerPixel, zu];

    Xdata_kriging = zeros(Ndrone, 2); 
    Ydata_kriging = zeros(Ndrone, 1);
    cnt = 0;
    ChannelModel.noise = 1;
    for idrone = 1:Ndrone
        DronePos = [droneXmat(idrone), droneYmat(idrone), DroneHeightSample];
        zu = 0;
        UserPos = [xu, yu, zu];
        gain = powersample([DronePos, UserPos], Maps, ChannelModel);
        cnt = cnt + 1;
        Xdata_kriging(cnt, :) = DronePos(1:2);
        Ydata_kriging(cnt, :) = gain;
    end
    
    d_variogram = variogram(Xdata_kriging, Ydata_kriging);
    [~, ~, ~, vstruct] = variogramfit(d_variogram.distance,d_variogram.val, ...
        [],[],[],'model', 'exponential', 'plotit', false);
    Gkri = kriging(vstruct, Xdata_kriging(:, 1), Xdata_kriging(:, 2), ...
        Ydata_kriging,meshgrid(1:lenX, 1:lenY)', meshgrid(1:lenY, 1:lenX));
    
    Gloc = powerMapReconBld(R,BldMapHat,FolMapHat,userPosAll(iuser, :), ...
        DroneHeightMap,Maps);
    Gknn = powerMapReconKNN(R,Ydata,userPosAll(iuser, :),DroneHeightMap,Maps);
    
    ChannelModel.noise = 0;
    UserPosArray = repmat([xu, yu, zu], length(DX(:)), 1);
    GainVec = powersample([DronePosArray UserPosArray], Maps, ChannelModel);% <--
    Gtru = reshape(GainVec, lenX, lenY);
    metrics_mse = [metrics_mse; mse(Gtru-Gkri) mse(Gtru-Gknn) mse(Gtru-Gloc)];
    metrics_mae = [metrics_mae; mean(abs(vec(Gtru-Gkri))) ...
        mean(abs(vec(Gtru-Gknn))) mean(abs(vec(Gtru-Gloc)))];
end
disp([mean(metrics_mae, 1) mean(metrics_mse, 1) Nue]);

close all
end
disp(datetime);
