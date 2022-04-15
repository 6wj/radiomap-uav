function Gr = powerMapReconKNN(R, Ydata, UE1PosMeter, DroneHeight, Maps)

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
        
        yx = 0;
        for n = 1:Nnear
            yx = yx + Ydata(I(n)) * exp( - d_vec(I(n)) ^ 2 / 2 / vard);
            % Zx = Zx + R.Z1(I(n), 1:K) * exp( - d_vec(I(n)) ^ 2 / 2 / vard);
            % Very heuristic reconstruction algorithm
            % Need to be optimized
        end
        yx = yx / (sum(exp( - d_vec(I(1:Nnear)) .^ 2 / 2 / vard)) + 1e-20);

%         %%% Modified
%         X1 = [DronePosMeter, UE1PosMeter];
%         Idx = knnsearch(R.X, X1);
%         yx = mean(Ydata(Idx));

        % 
        Gr(xi, yi) = yx;

    end
end