set(0,'defaultfigurecolor',[1 1 1])
clear
close all
addpath Sub

% Topology parameters
Nue = 100;       % Default 100
UE1 = [200 250]; %[200 250]; %[105, 120]; %[87, 123];           % meter

DroneHeightMap = 50;        % Drone height for air-to-ground power map
DroneHeightSample = 50;     % Drone height for learning (Learning stage)

noise = 4;
% Channel model 
C.A1 = - 22; % LOS (LTE TR36.814)
C.B1 = - 28;
C.S1 = noise;    % Model variance / shadowing in RMS (std)

C.A2 = - 28; % Obstructed LOS (foliage etc.)
C.B2 = - 24;
C.S2 = 3; 

C.A3 = - 36; % NLOS (LTE TR36.814)
C.B3 = - 22;
C.S3 = noise;

los_nlos_trans = 0;         % Size of transition region between segments
                            % This is to model the "edge effect", the
                            % signal should have a "smooth" transition from
                            % LOS to NLOS
ChannelModel.C = C;
ChannelModel.los_nlos_trans = los_nlos_trans;

% Urban topology data
% DATA = load('BldMap_trueHeight.mat');
% smallBldMap = ordfilt2(DATA.smallBldMap, 1, ones(3, 3));
% meterPerPixel = DATA.dsMeterPerPixel * 3;
% startp = 150;
% stopp = 450;
% BldMapZ = smallBldMap(startp:meterPerPixel:stopp, startp:meterPerPixel:stopp).';
% BldPosMat = (BldMapZ > 6);
DATA = load('BuildingMap_3m.mat');
meterPerPixel = DATA.dsMeterPerPixel;
BldMapZ = DATA.smallBldMap(1:end, 1:end).';
BldPosMat = DATA.BldAreaIndicator(1:end, 1:end).';

% % Topology data for test
% % Disable filter in advance
% BldMapZ = zeros(94, 94);
% BldMapZ(70, 40:60) = 30;
% BldMapZ(50, 40:60) = 35;
% BldMapZ(25, :) = 10;
% BldPosMat = zeros(94, 94);
% BldPosMat(70, 40:60) = 1;
% BldPosMat(50, 40:60) = 1;

% BldMapZ = zeros(94, 94);
% BldMapZ(70, 50) = 30;
% BldMapZ(50, 50) = 35;
% % BldMapZ(25, :) = 10;
% BldPosMat = zeros(94, 94);
% BldPosMat(70, 50) = 1;
% BldPosMat(50, 50) = 1;

[lenX, lenY] = size(BldMapZ);   % lenX = # of rows, corresponding to #
                                % of pixels in x-axis
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

%% DATA COLLECTION NONOISE ------------------------------------------------
ChannelModel.noise = 0;
UserPoseIndex0 = randperm(nUserPos);
UserPoseIndex = UserPosAll(UserPoseIndex0(1:Nue));
PosUE = zeros(Nue, 3);
RadioMap = zeros(Nue * lenX * lenY, 8);
for iuser = 1:Nue
    userIndex = UserPoseIndex(iuser);
    yu = floor((userIndex - 1) / lenX) + 1;
    xu = userIndex - (yu - 1) * lenX;
    zu = BldMapZ(xu, yu);
    UserPos = [xu, yu, zu];
    [DX, DY] = meshgrid((1:lenX), (1:lenY));
    DX = DX.';  DY = DY.';
    DronePosArray = [DX(:) DY(:) DroneHeightMap * ones(length(DX(:)), 1)];
    UserPosArray = repmat(UserPos, length(DX(:)), 1);
    GainVec = powersample([DronePosArray UserPosArray], Maps, ChannelModel);
    Gain = reshape(GainVec, lenX, lenY);
    DronePosMeter = [DronePosArray(:, 1:2) * meterPerPixel, DronePosArray(:, 3)];
    UserPosMeter = [UserPosArray(:, 1:2) * meterPerPixel, UserPosArray(:, 3)];
    PosUE(iuser, :) = [UserPos(1:2) * meterPerPixel, UserPos(3)];
    RadioMap((iuser - 1) * lenX * lenY + 1 : iuser * lenX * lenY, :) = ...
        [UserPosMeter, DronePosMeter, ...
        vecnorm(UserPosMeter - DronePosMeter, 2, 2), vec(Gain)];
end

save('radiomap_simulated100tx_3class.mat');
plot(log10(RadioMap(:, 7)), RadioMap(:, 8), '.');
