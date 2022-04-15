function Gr2 = powMapRecBldSig(R, BldMapHat, ~, UE1PosMeter, DroneHeight, Maps)
% powerMapRecon(R, UE1PosMeter, DroneHeight, lenX, lenY, meterPerPixel)

if nargin > 5
    [lenX, lenY] = size(Maps.BldMapZ);
    if size(BldMapHat, 1) < lenX || size(BldMapHat, 2) < lenY
        BldMapHat = upsamplematrix(BldMapHat, size(Maps.BldPosMat));
    end
    meterPerPixel = Maps.meterPerPixel;

    BldMapHat(Maps.BldPosMat < 1 & Maps.FolPosMat < 1) = UE1PosMeter(3);
% Note: The actual ground level may not be 0 meter. In practice, we do not
% assume the knowledge of altitude of the ground level; instead, we only
% need to compare the altitude at every foliage and building location
% (assumed known). However, for coding coveniece, we compare the altitude
% at every location (including non-foliage, non-building locations), and
% hence we should shift the ground level (only of non-foliage, non-building 
% locations). 
else
    [lenX, lenY] = size(BldMapHat);
    meterPerPixel = 1;
end

UE1_pixel = max(round(UE1PosMeter(1:2) / meterPerPixel), 1);
UE1Pos = [UE1_pixel, UE1PosMeter(3)];
Gr2 = zeros(lenX, lenY);
for xi = 1:lenX
    for yi = 1:lenY
        DronePos = [xi, yi, DroneHeight];   % [pixel, pixel, meter]
        DronePosMeter = [DronePos(1:2) * meterPerPixel, DronePos(3)];
        d = norm(UE1PosMeter - DronePosMeter, 2);

        noBld = 1;
        [covBlds, covBldZs] = covBldZ(DronePos, UE1Pos, lenX);
        for i = 1:length(covBlds)
            j = covBlds(i);
            yb = floor((j - 1) / lenX) + 1;
            xb = j - (yb - 1) * lenX;
            zj = covBldZs(i);
            noBld = noBld * sigmoid(zj - BldMapHat(xb, yb));
        end

        K = length(R.Alpha);
        glos = R.Alpha(1) * log10(d) + R.Beta(1);
        gnlos = R.Alpha(K) * log10(d) + R.Beta(K);
        Gr2(xi, yi) = glos * noBld + gnlos * (1 - noBld);

    end
end