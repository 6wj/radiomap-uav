function [Hs, S, Hs1, S1] = covPosZ(R, Maps, lenX, lenY)
% Measurements pass obstacles
%
% INPUT
%   R           R from EM algorithm
%   Map         Reconstructed city map 
%
% OUTPUT
%   Hs          Index set of buildings under the line segment
%   S

N = size(R.X, 1);

% - Building
BldPos = find(Maps.BldPosMat > 0);
nBldPos = length(BldPos);
Hs = cell(N, 2);        % cell array for index and Z-value
S = cell(nBldPos, 1);   % the j-th cell contains the set of measurements that invole the j-th building
for idata = 1:N
    DronePosMeter = R.X(idata, 1:3);
    DronePosPixel = [floor(DronePosMeter(1:2) / Maps.meterPerPixel), DronePosMeter(3)];
    UserPosMeter = R.X(idata, 4:6);
    UserPosPixel = [floor(UserPosMeter(1:2) / Maps.meterPerPixel), UserPosMeter(3)];
    
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
    end
end

% - Foliage
FoliagePos = find(Maps.FolPosMat > 0);
nFoliagePos = length(FoliagePos);

Hs1 = cell(N, 2);        % cell array for index and Z-value
S1 = cell(nFoliagePos, 1);   % the j-th cell contains the set of measurements that invole the j-th building
for idata = 1:N
    DronePosMeter = R.X(idata, 1:3);
    DronePosPixel = [floor(DronePosMeter(1:2) / Maps.meterPerPixel), DronePosMeter(3)];
    UserPosMeter = R.X(idata, 4:6);
    UserPosPixel = [floor(UserPosMeter(1:2) / Maps.meterPerPixel), UserPosMeter(3)];
    
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
    end
end

end