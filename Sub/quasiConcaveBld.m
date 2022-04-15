function [BldMapHat, FolMapHat] = quasiConcaveBld(R, Maps, ...
    lenX, lenY, DroneHeightMap, H0, Hs, S, Hs1, S1)
% Quasi-concave height estimation 
% 
% INPUT
%   R           Output from the EM algorithm 
%   Maps        A structure that contains: 
%       Maps.BldPosMat       An indicator matrix for the positions of buildings
%       Maps.FolPosMat       An indicator matrix for the positions of trees
%       Maps.meterPerPixel   A scalar for how many meters for a matrix
%                            element

%global metrics;

% Connected-component labeling for building blobs classification
divideSize = round(4 * Maps.meterPerPixel);
blobsNum = 0;
blobsLabel = zeros(lenX, lenY);
for j = 1:divideSize:lenY
    for i = 1:divideSize:lenX
% lenX = # of rows, corresponding to # of pixels in x-axis
% lenY = # of columns, corresponding to # of pixels in y-axis
        [blobsLabelPart, blobsNumPart] = bwlabel(...
            Maps.BldPosMat(i:min(i+divideSize-1,lenX),j:min(j+divideSize-1,lenY)));
        blobsLabel(i:min(i+divideSize-1,lenX),j:min(j+divideSize-1,lenY))...
            = blobsLabelPart + blobsNum;
        blobsNum = blobsNum + blobsNumPart;
    end
end
blobsLabel = blobsLabel .* Maps.BldPosMat;
blobsLabel = blobsLabel(Maps.BldPosMat > 0);
nBldBlob = cell(blobsNum, 2);
mapsNum = 1:length(blobsLabel);
for i = 1:blobsNum
    nBldBlob{i, 1} = mapsNum(blobsLabel == i);
    nBldBlob{i, 2} = mapsNum(blobsLabel > 0 & blobsLabel ~= i);
end

% - Building height
maxBldHeight = DroneHeightMap;  % meter
BldPos = find(Maps.BldPosMat > 0);

Ps = R.Z1; % data for los/nlos probability

%
show_figures = 1;
B0 = H0(Maps.BldPosMat > 0) * maxBldHeight;
B0(B0 >= maxBldHeight) = maxBldHeight;
Nbld = length(B0(:));
eps_h = 0.001;
MAXLOOP = 7;      % maximum number of iterations
H_array = zeros(MAXLOOP, Nbld);
cnt = 1;
H_array(cnt, :) = B0;

% Xtr = [reshape((ones(lenX, 1)*(1:lenY))', [lenX * lenY 1]) ...
%     reshape((ones(lenX, 1)*(1:lenY)), [lenX * lenY 1])];
% Ytr = reshape(BldMapHat, [lenX * lenY 1]);
% predFunc = localPolyRegressionCV(Xtr, H', 2, 1, kernelParams);

while cnt <= 3 || cnt < MAXLOOP && ...
        norm(H_array(cnt, :) - H_array(cnt - 1, :)) > eps_h * Nbld
    if show_figures
        fprintf('Building height estimation iteration: %2d', cnt);
    end
    
    if cnt <= 2
        dsfactor = 1;
    elseif cnt <= 3
        dsfactor = 1;
    else
        dsfactor = 1;
    end
    
    %metrics = zeros(1, length(Ps));
    H = zeros(1, Nbld);
    H_array_cnt = H_array(cnt, :);
    parfor k = 1:Nbld
        % Estimate the k-th height
        H(k) = mlHeight(k, H_array_cnt, Hs, Ps, S, dsfactor, ...
            maxBldHeight, nBldBlob, blobsLabel);
    end
    if show_figures
        fprintf('\n');
    end
    %metrics_mat = [metrics_mat sum(metrics)];
    cnt = cnt + 1;
    H_array(cnt, :) = H;
    
end
BldMapHat = zeros(size(Maps.BldPosMat));
BldMapHat(BldPos) = H_array(cnt, :);
BldMapHat = medfilt2(BldMapHat .* (BldMapHat <= maxBldHeight)) .* Maps.BldPosMat;

% - Foliage height
maxFolHeight = 25;  % meter
FoliagePos = find(Maps.FolPosMat > 0);

%
show_figures = 1;
F0 = H0(Maps.FolPosMat > 0) * maxFolHeight;
F0(F0 >= maxFolHeight) = maxFolHeight;
Nbld = length(F0(:));
eps_h = 0.001;
MAXLOOP = 8;      % maximum number of iterations
H_array1 = zeros(MAXLOOP, Nbld);
cnt = 1;
H_array1(cnt, :) = F0;

while cnt <= 3 || cnt < MAXLOOP && ...
        norm(H_array1(cnt, :) - H_array1(cnt - 1, :)) > eps_h * Nbld
    if show_figures
        fprintf('Foliage height estimation iteration:  %2d', cnt);
    end
    
    H = zeros(1, Nbld);
    H_array_cnt = H_array1(cnt, :);
    parfor k = 1:Nbld
        % Estimate the k-th height
        H(k) = mlHeight(k, H_array_cnt, Hs1, Ps, S1, 1, maxFolHeight);
    end
    if show_figures
        fprintf('\n');
    end
    cnt = cnt + 1;
    H_array1(cnt, :) = H;
    
end
FolMapHat = zeros(size(Maps.FolPosMat));
FolMapHat(FoliagePos) = H_array1(cnt, :);
FolMapHat = FolMapHat .* (FolMapHat <= maxFolHeight);

end