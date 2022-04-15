function [smallMaps, cols, covB] = downsampleMaps(Maps, dsfactor, R)

smallMaps = Maps;
smallMaps.BldMapZ = dsmap(Maps.BldMapZ, dsfactor);
smallMaps.BldPosMat = dsmap(Maps.BldPosMat, dsfactor);
smallMaps.FolPosMat = dsmap(Maps.FolPosMat, dsfactor);
smallMaps.meterPerPixel = Maps.meterPerPixel * dsfactor;

if ~isfield(Maps, 'neighbourMat')
    return
end
[lenX, lenY] = size(smallMaps.BldPosMat);
meterPerPixel = smallMaps.meterPerPixel;
droneHeightMap = Maps.droneHeightMap;
neighbourMat = Maps.neighbourMat;
obstacles = (1:lenX*lenY)';
nObst = length(obstacles);
nMeas = length(R.X);
cols = cell(nObst, 1);
covB = cell(nMeas, size(neighbourMat, 1) * 2);
for nr = 1:size(neighbourMat, 1)
    for id = 1:nMeas
        droneMeter = R.X(id, 1:3);
        droneMeter = [droneMeter(1:2)+neighbourMat(nr, :) droneMeter(3)];
        dronePixel = [floor(droneMeter(1:2)/meterPerPixel)+1, droneMeter(3)];
        dronePixel = max(dronePixel, 1);
        dronePixel = min(dronePixel, [lenX lenY droneMeter(3)]);
        userMeter = R.X(id, 4:6);
        userPixel = [floor(userMeter(1:2)/meterPerPixel)+1, userMeter(3)];
        [covBlds, covBldZs] = covBldZ(dronePixel, userPixel, lenX);
        covBlds = covBlds(covBldZs > 0 & covBldZs <= droneHeightMap);
        covBldZs = covBldZs(covBldZs > 0 & covBldZs <= droneHeightMap);
        covBlds = covBlds(covBlds > 0 & covBlds <= lenX * lenY);
        covBldZs = covBldZs(covBlds > 0 & covBlds <= lenX * lenY);
        covB{id, nr * 2 - 1} = covBlds;
        covB{id, nr * 2} = covBldZs;
        for i = 1:length(covBlds)
            j = covBlds(i);
            ib = find(obstacles == j, 1, 'first');
            if ~isempty(ib)
                if ~isempty(cols{ib})
                    cols{ib} = unique([cols{ib} id]);
                else
                    cols{ib} = id;
                end
            end
        end
    end
end

end


function smallA = dsmap(A, dsfactor2)

% dsfactor2 = 4; % 18.7617; % 18.7617 -> 3 meter precision
[lenX, lenY] = size(A);

smallLenX = ceil(lenX / dsfactor2);
smallLenY = ceil(lenY / dsfactor2);
smallA = zeros(smallLenX, smallLenY);

for i = 1:smallLenX
    for j = 1:smallLenY
        xij = (i - 1) * dsfactor2 + 1;
        yij = (j - 1) * dsfactor2 + 1;
        x0 = max(1, min(lenX, round(xij)));
        y0 = max(1, min(lenY, round(yij)));
        x1 = max(1, min(lenX, round(xij + dsfactor2)));
        y1 = max(1, min(lenY, round(yij + dsfactor2)));
        subblock = A(x0:x1, y0:y1);
        smallA(i, j) = mean(subblock(:));
    end
end
end