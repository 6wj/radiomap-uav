function [Gain, noBld, noFoliage] = powersample(PosArray, Maps, ChannelModel)
%
% PosArray = [DronePosArray, UserPosArray]
%          = [Pixel, Pixel, meter, pixel, pixel, meter];
%
SMOOTH_TRANSITION= 1;

BldMapZ = Maps.BldMapZ;
BldPosMat = Maps.BldPosMat;
FolPosMat = Maps.FolPosMat;
meterPerPixel = Maps.meterPerPixel;

C = ChannelModel.C;
noise = ChannelModel.noise;
los_nlos_trans = ChannelModel.los_nlos_trans;

if los_nlos_trans < 1e-9
    los_nlos_trans = 1e-9;
end

[lenX, ~] = size(BldMapZ);
[N, D] = size(PosArray);
dim = round(D / 2);

Gain = zeros(1, N);
for ipos = 1:N
    DronePos = PosArray(ipos, 1:dim);
    UserPos = PosArray(ipos, dim + 1:end);
    
    noBld = 1;
    noFoliage = 1;
    dHeightBld = -1;
    dHeightFol = -1;
    [covBlds, covBldZs] = covBldZ(DronePos, UserPos, lenX);
    for i = 1:length(covBlds)
        j = covBlds(i);
        yb = floor((j - 1) / lenX) + 1;
        xb = j - (yb - 1) * lenX;
        zj = covBldZs(i);
        if zj < BldMapZ(xb, yb)  
            if BldPosMat(xb, yb) == 1
                noBld = 0;
                if BldMapZ(xb, yb) - zj > dHeightBld
                    dHeightBld = BldMapZ(xb, yb) - zj;
                end
            end 
            if FolPosMat(xb, yb) == 1
                noFoliage = 0;
                if BldMapZ(xb, yb) - zj > dHeightFol
                    dHeightFol = BldMapZ(xb, yb) - zj;
                end
            end

        end
    end

    DronePosMeter = [DronePos(1:2) * meterPerPixel, DronePos(3)];
    UserPosMeter = [UserPos(1:2) * meterPerPixel, UserPos(3)];
    d = norm(UserPosMeter - DronePosMeter, 2);

    if noise > 0
        glos = C.A1 * log10(d) + C.B1 + randn * sqrt(C.S1);
        golos = C.A2 * log10(d) + C.B2 + randn * sqrt(C.S2);
        gnlos = C.A3 * log10(d) + C.B3 + randn * sqrt(C.S3);
    else
        glos = C.A1 * log10(d) + C.B1;
        golos = C.A2 * log10(d) + C.B2;
        gnlos = C.A3 * log10(d) + C.B3;
    end
    if noBld && noFoliage
        % LOS
        Gain(ipos) = glos;    % CAUTION about the matrix form
    elseif noBld && ~ noFoliage
        % Obstructed LOS
        if SMOOTH_TRANSITION
            % Smooth transition
            if dHeightFol < los_nlos_trans
                Gain(ipos) = (los_nlos_trans - dHeightFol) / los_nlos_trans * glos ...
                                + dHeightFol / los_nlos_trans * golos;
            else
                Gain(ipos) = golos;
            end
        else
            % Hard transition
            Gain(ipos) = golos;
        end
    else
        % NLOS
        if SMOOTH_TRANSITION
            % Smooth transition
            if dHeightBld < los_nlos_trans
                Gain(ipos) = (los_nlos_trans - dHeightBld) / los_nlos_trans * golos ...
                                + dHeightBld / los_nlos_trans * gnlos;
            else
                Gain(ipos) = gnlos;
            end
        else
            % Hard transition
            Gain(ipos) = gnlos;
        end

    end

end