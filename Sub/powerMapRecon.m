function Gr = powerMapRecon(R, UE1PosMeter, DroneHeight, Maps)
% powerMapRecon(R, UE1PosMeter, DroneHeight, lenX, lenY, meterPerPixel)

[lenX, lenY] = size(Maps.BldMapZ);
meterPerPixel = Maps.meterPerPixel;

K = length(R.Alpha);
Ndat = size(R.X, 1);
% UE1_pixel = [1, 1] + floor(UE1 / meterPerPixel);
% UE1Pos = [UE1_pixel, BldMapZ(UE1_pixel(2), UE1_pixel(1))];          % [pixel, pixel, meter]
% UE1PosMeter = [UE1, UE1Pos(3)];
% DroneHeight = 50;           % meter
Gr = zeros(lenX, lenY);
for xi = 1:lenX
    for yi = 1:lenY
        DronePos = [xi, yi, DroneHeight];   % [pixel, pixel, meter]
        DronePosMeter = [DronePos(1:2) * meterPerPixel, DronePos(3)];
        d = norm(UE1PosMeter - DronePosMeter, 2);

        %
        Nnear = 5;
        X1 = [DronePosMeter, UE1PosMeter];
        Zx = zeros(1, K);
        d_vec = zeros(Ndat, 1);
        for l = 1:Ndat
            d_vec(l) = norm(X1 - R.X(l, :));
        end
        [~, I] = sort(d_vec, 'ascend');
        % vard = var(d_vec(I(1:Nnear)));
        vard = mean(d_vec(I(1:Nnear)) .^ 2);
        for n = 1:Nnear
            Zx = Zx + R.Z1(I(n), 1:K) * exp( - d_vec(I(n)) ^ 2 / 2 / vard);
            % Very heuristic reconstruction algorithm
            % Need to be optimized
        end
        Zx = Zx / (1e-10 + sum(Zx));

        yx = 0;
        for k = 1:K
            if Zx(k) > 0
                yx = yx + Zx(k) * (R.Alpha(k) * log10(d) + R.Beta(k));
            end
        end
        % 
        Gr(xi, yi) = yx;

    end
end