function Gr2 = powerMapReconKBld(R, H, UE1PosMeter, DroneHeight, Maps)

meterPerPixel = Maps.meterPerPixel;
UE1_pixel = max(round(UE1PosMeter(1:2) / meterPerPixel), 1);
UE1Pos = [UE1_pixel, UE1PosMeter(3)];
[lenX, lenY] = size(Maps.BldMapZ);
Gr2 = zeros(lenX, lenY);
K = length(R.Alpha) - 1;
for xi = 1:lenX
    for yi = 1:lenY
        DronePos = [xi, yi, DroneHeight];   % [pixel, pixel, meter]
        DronePosMeter = [DronePos(1:2) * meterPerPixel, DronePos(3)];
        d = norm(UE1PosMeter - DronePosMeter, 2);

        [covBlds, covBldZs] = covBldZ(DronePos, UE1Pos, lenX);
        for k = 0:K
            if sum(vec(H(k+1:K, covBlds) > ...
                    repmat(covBldZs, K-k, 1))) == 0 && ...
                    (k == 0 || ...
                    sum(H(k, covBlds) > covBldZs) > 0)
                Gr2(xi, yi) = (R.Alpha(k+1) * log10(d) + R.Beta(k+1));
            end
        end
    end
end