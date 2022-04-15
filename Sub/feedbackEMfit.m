function R = feedbackEMfit(R, X, Y, Map, BldPosMat, FolPosMat, meterPerPixel)
% Segmented regression for the channel model based on maximul likelihood
% (ML) estimation and expectation-maximizaiton (EM) algorithm
%
% INPUT
%   R0          R from last regression EM algorithm
%   Map         Reconstructed city map 
%
% OUTPUT
%   R.Alpha
%   R.Beta
%   R.Sigma2
%   R.Pi
%   R.Xc

maxSigma1 = inf;
N = length(X);

Xdata = R.X;
[lenX, ~] = size(Map);
for idata=1:N
    DronePosMeter = Xdata(idata, 1:3);
    DronePos = [floor(DronePosMeter(1:2) / meterPerPixel), DronePosMeter(3)];
    UserPosMeter = Xdata(idata, 4:6);
    UserPos = [floor(UserPosMeter(1:2) / meterPerPixel), UserPosMeter(3)];
    noBld = 1;
    noFoliage = 1;
    [covBlds, covBldZs] = covBldZ(DronePos, UserPos, lenX);
    for i = 1:length(covBlds)
        j = covBlds(i);
        yb = floor((j - 1) / lenX) + 1;
        xb = j - (yb - 1) * lenX;
        zj = covBldZs(i);
        if zj < Map(xb, yb)  
            if BldPosMat(xb, yb) == 1
                noBld = 0;
            end 
            if FolPosMat(xb, yb) == 1
                noFoliage = 0;
            end

        end
    end
    if noBld && noFoliage
        % LOS
        R.Z1(idata, :) = 1;
    elseif noBld && ~ noFoliage
        % Obstructed LOS
        R.Z1(idata, :) = 0;
    else
        % NLOS
        R.Z1(idata, :) = 0;
    end
end

I1 = R.Z1 > 0.5; I2 = R.Z1 < 0.5;

maxloop = 100; cnt = 0; f = 1; f0 = 0; F = zeros(1, maxloop);
% EM loop begins --
while abs(f - f0) > 1e-9 && cnt < maxloop
    cnt = cnt + 1; f0 = f;
    
    % Regression step
    P1 = polyfit(X(I1), Y(I1), 1);
    P2 = polyfit(X(I2), Y(I2), 1);
    Z1 = Y - (P1(1) * X + P1(2));     % Distance for all data points to model 1
    Z2 = Y - (P2(1) * X + P2(2));     % Distance for all data points to model 2
    S1 = var(Z1(I1)); S2 = var(Z2(I2)); % New variance
    S1 = min(maxSigma1^2, S1);          % Our heuristic 1: LOS variance upper bound
                                        % however, this does not work
                                        % (hence disabled)

    % Classification
    D1 = Z1.^2 / S1; D2 = Z2.^2 / S2;       % Normalized squared distances
    L1 = exp(- D1 / 2); L2 = exp(- D2 / 2); % Likelihood 
    C1 = L1 ./ (L1 + L2);                   % Normalized likelihood to model 1 (LOS)
                                            % Note, the soft label formula is 
                                            % different from that in the
                                            % paper

    gamma = 0.5; 
    C1(Z1 > 0) = 1;                     % Heuristic 2: Modify the Bayesian 
    I1 = C1 > gamma;                    % probability such that NLOS cannot
    I2 = ~ (C1 > gamma);                % be stronger than LOS
    
    f = sum(C1 .* D1 + (1 - C1) .* D2) / N; % Objective: clustered distance
    F(cnt) = f; 
    
end

Ap = [P1(1), P2(1)]; Bp = [P1(2), P2(2)];

R.Alpha = Ap;
R.Beta = Bp;
R.Sigma = sqrt([S1, S2]);
R.Z1 = C1;

end