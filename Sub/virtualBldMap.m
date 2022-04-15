function [BldMapHat, FolMapHat, w] = virtualBldMap(R, Maps, H0)
% Building/Foliage height estimation 
% [BldMapHat, FolMapHat] = virtualBldMap(R, Maps, Map0)
% 
% INPUT
%   R           Output from the EM algorithm 
%   Maps        A structure that contains: 
%       Maps.BldPosMat       An indicator matrix for the positions of buildings
%       Maps.FolPosMat       An indicator matrix for the positions of trees
%       Maps.meterPerPixel   A scalar for how many meters for a matrix
%                            element
%   H0          Initial building height values
 
klos = 1;
kolos = 2;
knlos = 3;
length_R_Z1 = size(R.Z1, 2);
N = size(R.X, 1);   % number of data samples
meterPerPixel = Maps.meterPerPixel;
lenX = size(Maps.BldPosMat, 1);

% - Building height
maxBldHeight = 50;
BldPos = find(Maps.BldPosMat > 0);
nBldPos = length(BldPos);
if nargin < 3
    B0 = ones(1, nBldPos) * maxBldHeight;
else 
    B0 = H0(Maps.BldPosMat > 0);
    B0(B0 > maxBldHeight) = maxBldHeight;
end

Ps = zeros(N, 1);       % data for los probability
Hs = cell(N, 2);        % cell array for index and Z-value
S = cell(nBldPos, 1);   % the j-th cell contains the set of measurements that invole the j-th building
for idata = 1:N
    DronePosMeter = R.X(idata, 1:3);
    DronePosPixel = [1 1 0] + [floor(DronePosMeter(1:2) / meterPerPixel), DronePosMeter(3)];
    UserPosMeter = R.X(idata, 4:6);
    UserPosPixel = [1 1 0] + [floor(UserPosMeter(1:2) / meterPerPixel), UserPosMeter(3)];
    
    [covBlds, covBldZs] = covBldZ(DronePosPixel, UserPosPixel, lenX);
    BldIds = zeros(1, length(covBlds));
    BldZval = zeros(1, length(covBlds));
    cnt = 0;
    for i = 1:length(covBlds)
        j = covBlds(i);
        ib = find(BldPos == j, 1, 'first');
        if ~isempty(ib)
            cnt = cnt + 1;
            BldIds(cnt) = ib;
            BldZval(cnt) = covBldZs(i);
            if ~isempty(S{ib})
                S{ib} = [S{ib} idata];
            else
                S{ib} = idata;
            end
        end
    end
    if cnt > 0
        Hs{idata, 1} = BldIds(1:cnt);
        Hs{idata, 2} = BldZval(1:cnt);
        if length_R_Z1 == 5
            Ps(idata) = R.Z1(idata, klos) + R.Z1(idata, kolos);
        elseif length_R_Z1 ~= 2
            Ps(idata) = R.Z1(idata, klos);
        end
    end
end
if length_R_Z1 == 2
    Ps = R.Z1;
end

%
Hhat = bldHeight3(Hs, Ps, S, maxBldHeight, B0, 0);
BldMapHat = zeros(size(Maps.BldPosMat));
BldMapHat(BldPos) = Hhat;

% - Foliage height
maxFolHeight = 25;  % meter
FoliagePos = find(Maps.FolPosMat > 0);
nFoliagePos = length(FoliagePos);

if nargin < 3
    F0 = ones(1, nFoliagePos) * maxFolHeight;
else 
    F0 = H0(Maps.FolPosMat > 0);
    F0(F0 > maxFolHeight) = maxFolHeight;
end

Ps1 = zeros(N, 1);          % data for los probability
Hs1 = cell(N, 2);        % cell array for index and Z-value
S1 = cell(nFoliagePos, 1);   % the j-th cell contains the set of measurements that invole the j-th building
for idata = 1:N
    DronePosMeter = R.X(idata, 1:3);
    DronePosPixel = [1 1 0] + [floor(DronePosMeter(1:2) / meterPerPixel), DronePosMeter(3)];
    UserPosMeter = R.X(idata, 4:6);
    UserPosPixel = [1 1 0] + [floor(UserPosMeter(1:2) / meterPerPixel), UserPosMeter(3)];
    
    [covFols, covFolZs] = covBldZ(DronePosPixel, UserPosPixel, lenX);
    FolIds = zeros(1, length(covFols));
    FolZval = zeros(1, length(covFols));
    cnt = 0;
    for i = 1:length(covFols)
        j = covFols(i);
        ib = find(FoliagePos == j, 1, 'first');
        if ~isempty(ib)
            cnt = cnt + 1;
            FolIds(cnt) = ib;
            FolZval(cnt) = covFolZs(i);
            if ~isempty(S1{ib})
                S1{ib} = [S1{ib} idata];
            else
                S1{ib} = idata;
            end
        end
    end
    if cnt > 0
        Hs1{idata, 1} = FolIds(1:cnt);
        Hs1{idata, 2} = FolZval(1:cnt);
        Ps1(idata) = R.Z1(idata, klos);
    end
end

%
maxFolHeight = 25;
Folhat = bldHeight3(Hs1, Ps1, S1, maxFolHeight, F0, 0);
FolMapHat = zeros(size(Maps.FolPosMat));
FolMapHat(FoliagePos) = Folhat;

% Generate w, LoS/NLoS labels from heights
w = true(N, 1);
for i = 1:N
    w(i) = isempty(Hs{i, 1}) | sum(vec(Hhat(Hs{i, 1})) > vec(Hs{i, 2})) == 0;
end