function Gr2 = powerMapReconBld(R, BldMapHat, FolMapHat, UE1PosMeter, DroneHeight, Maps, PosMatOn)
% powerMapRecon(R, UE1PosMeter, DroneHeight, lenX, lenY, meterPerPixel)

if nargin > 5
    [lenX, lenY] = size(Maps.BldMapZ);
    if size(BldMapHat, 1) < lenX || size(BldMapHat, 2) < lenY
        BldMapHat = upsamplematrix(BldMapHat, size(Maps.BldPosMat));
        FolMapHat = upsamplematrix(FolMapHat, size(Maps.BldPosMat));
    end
    meterPerPixel = Maps.meterPerPixel;

    % BldMapHat(Maps.BldPosMat < 1 & Maps.FolPosMat < 1) = Maps.BldMapZ(Maps.BldPosMat < 1 & Maps.FolPosMat < 1);
    % FolMapHat(Maps.BldPosMat < 1 & Maps.FolPosMat < 1) = Maps.BldMapZ(Maps.BldPosMat < 1 & Maps.FolPosMat < 1);
    BldMapHat(Maps.BldPosMat < 1 & Maps.FolPosMat < 1) = UE1PosMeter(3);
    FolMapHat(Maps.BldPosMat < 1 & Maps.FolPosMat < 1) = UE1PosMeter(3);
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
        noFoliage = 1;
        [covBlds, covBldZs] = covBldZ(DronePos, UE1Pos, lenX);
        for i = 1:length(covBlds)
            j = covBlds(i);
            yb = floor((j - 1) / lenX) + 1;
            xb = j - (yb - 1) * lenX;
            zj = covBldZs(i);
            if zj < BldMapHat(xb, yb)  
                noBld = 0;
            end
            if zj < FolMapHat(xb, yb)
                noFoliage = 0;
            end
        end

        if exist('PosMatOn', 'var')
            glos = true;
            golos = true;
            gnlos = false;
        else
            K = length(R.Alpha);
            glos = R.Alpha(1) * log10(d) + R.Beta(1);
            golos = R.Alpha(2) * log10(d) + R.Beta(2);
            gnlos = R.Alpha(K) * log10(d) + R.Beta(K);
        end
        if noBld && noFoliage
            Gr2(xi, yi) = glos;    
        elseif noBld && ~ noFoliage
            Gr2(xi, yi) = golos;
        else
            Gr2(xi, yi) = gnlos; 
        end

    end
end