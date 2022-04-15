function R = segregress2K(G, Y, S, R, K, maxVar)
% Version 2.1 (Change the initialization method)
% R = segregress(X, Y, K)
%
% Segmented regression for the channel model based on maximul likelihood
% (ML) estimation and expectation-maximizaiton (EM) algorithm
%
% INPUT
%   X           N * D matrix, for N data samples; each row contains a
%               sample with d dimension. 
%   Y           N * 1 vector, for N data samples. 
%   K           The number of segments/clusters/regions
%
% OUTPUT
%   R.Alpha
%   R.Beta
%   R.Sigma2
%   R.Pi
%   R.Xc

if nargin < 6
    maxVar = single(intmax);
end

N = length(G);
G2 = G .^ 2;

% Initialization
Z1 = zeros(N, K * 2 - 1);
for i = 1:K * 2 - 1
    Z1(:, i) = (S == i);
end
Gtest = max(G);
Ytest = R.Alpha * Gtest + R.Beta;
[~, I] = sort(Ytest, 'descend');
Alpha = R.Alpha(I);
Beta = R.Beta(I);
Sigma2 = R.Sigma2(I);
Z1 = Z1(:, I);
Pi = R.Pi(I);

% EM algorithm ------------------------------------------------------------
Z0 = zeros(N, K + K - 1);   % K propagation segments + K - 1 transit segments
Z1 = [Z1 zeros(size(Z1, 1), K - 1)];
Pi = [Pi * 0.9, ones(1, K - 1) * 0.1 / (K - 1)];
% Sigma2(Sigma2 > maxVar) = maxVar;
Sigma2 = [Sigma2, maxVar * ones(1, K - 1)];
Alpha = [Alpha, zeros(1, K - 1)];
Beta = [Beta, zeros(1, K - 1)];
eps = 1e-3;
i0 = i;
i = 0;

MAXLOOP = 2;
% Q = zeros(1, MAXLOOP);
while norm(Z0 - Z1, 'fro') / norm(Z1, 'fro') > eps && i < MAXLOOP
    i = i + 1;
    
    Z0 = Z1;
    % Handle the variance
    for k = 1:K
        if Sigma2(k) > maxVar

            [~, Im] = max(Z1, [], k);
            Ik = find(Im == k);
            
            Dev_k = zeros(1, length(Ik));
            for cnt = 1:length(Ik)
                if k >= 2 && Y(Ik(cnt)) > Alpha(k) * G(Ik(cnt)) + Beta(k)
                      
                      Dev_k(cnt) = (Y(Ik(cnt)) - (Alpha(k) * G(Ik(cnt)) + Beta(k))) ...
                                / (Alpha(k - 1) * G(Ik(cnt)) + Beta(k - 1) ...
                                        - (Alpha(k) * G(Ik(cnt)) + Beta(k)));
                            
                elseif k < K && Y(Ik(cnt)) < Alpha(k) * G(Ik(cnt)) + Beta(k)
                    
                    Dev_k(cnt) = (Alpha(k) * G(Ik(cnt)) + Beta(k) - Y(Ik(cnt))) ...
                               / (Alpha(k) * G(Ik(cnt)) + Beta(k) ...
                                        - (Alpha(k + 1) * G(Ik(cnt)) + Beta(k + 1)));
                end                
            end
            
            [~, Iksort] = sort(Dev_k, 'descend');
            Izsort = Ik(Iksort);
            
            A = zeros(2);
            b = zeros(2, 1);
            A(1, 1) = sum(Z1(Izsort, k) .* G2(Izsort));
            A(1, 2) = sum(Z1(Izsort, k) .* G(Izsort));
            A(2, 1) = A(1, 2);
            A(2, 2) = sum(Z1(Izsort, k));
            b(1) = sum(Z1(Izsort, k) .* G(Izsort) .* Y(Izsort));
            b(2) = sum(Z1(Izsort, k) .* Y(Izsort));

            cnt = 0;
            while Sigma2(k) > maxVar && cnt < length(Izsort)
                cnt = cnt + 1;
                if k >= 2 && Y(Izsort(cnt)) > Alpha(k) * G(Izsort(cnt)) + Beta(k) ...
                          && Y(Izsort(cnt)) < Alpha(k - 1) * G(Izsort(cnt)) + Beta(k - 1)
                      
                    A(1, 1) = A(1, 1) - Z1(Izsort(cnt), k) * G2(Izsort(cnt)); 
                    A(1, 2) = A(1, 2) - Z1(Izsort(cnt), k) * G(Izsort(cnt));
                    A(2, 1) = A(1, 2);
                    A(2, 2) = A(2, 2) - Z1(Izsort(cnt), k);
                    b(1) = b(1) - Z1(Izsort(cnt), k) * G(Izsort(cnt)) * Y(Izsort(cnt));
                    b(2) = b(2) - Z1(Izsort(cnt), k) * Y(Izsort(cnt));
                    
                    Z1(Izsort(cnt), K + k - 1) = Z1(Izsort(cnt), k);
                    Z1(Izsort(cnt), k) = 0;
                    
                elseif k < K && Y(Izsort(cnt)) < Alpha(k) * G(Izsort(cnt)) + Beta(k)
                    
                    A(1, 1) = A(1, 1) - Z1(Izsort(cnt), k) * G2(Izsort(cnt)); 
                    A(1, 2) = A(1, 2) - Z1(Izsort(cnt), k) * G(Izsort(cnt));
                    A(2, 1) = A(1, 2);
                    A(2, 2) = A(2, 2) - Z1(Izsort(cnt), k);
                    b(1) = b(1) - Z1(Izsort(cnt), k) * G(Izsort(cnt)) * Y(Izsort(cnt));
                    b(2) = b(2) - Z1(Izsort(cnt), k) * Y(Izsort(cnt));
                    
                    Z1(Izsort(cnt), K + k) = Z1(Izsort(cnt), k);
                    Z1(Izsort(cnt), k) = 0;
                    
                end
                
                if Z1(Izsort(cnt), k) < 1e-9 &&...
                        ~isnan(trace(A)) && ~isinf(trace(A)) && rank(A) == 2

                    alpha_beta = A \ b;
                    
                    Alpha(k) = alpha_beta(1);
                    Beta(k) = alpha_beta(2);
                    Sigma2(k) = sum( Z1(Ik, k) .* (Y(Ik) - Alpha(k) * G(Ik) - Beta(k)).^2 ) ...
                                / sum( Z1(Ik, k)); 

%                     % DEBUG        
%                     R.Alpha = Alpha;
%                     R.Beta = Beta;
%                     R.Sigma2 = Sigma2;
%                     R.Pi = Pi;
%                     R.Q = Q(1:i);
%                     R.Z1 = Z1;
%                     R.X = X;
%                     debug_segregress(X, Y, R);
                end
            end
        end
    end
    
    % E-step: Update the soft label (i.e., the probably Z(l, k) that the l-th
    % sample is assigned to the k-th cluster
    % Z0 = Z1;
    Z1 = zeros(N, K + K - 1);
    Pxy = zeros(N, K + K - 1);      % The joint density p_k(x, y) for the k-th cluster eval @ l-th data sample (x, y)
    for l = 1:N
        for k = 1:K
            Pxy(l, k) = exp( - (Y(l) - Alpha(k) * G(l) - Beta(k))^2 / 2 / Sigma2(k)) ...
                        / (sqrt(2 * pi * Sigma2(k)));
        end
        for k = 1:K - 1 % For transit segments
            if i == 1
                alow = Alpha(k + 1) * G(l) + Beta(k + 1) + sqrt(Sigma2(k + 1));
                ahigh = Alpha(k) * G(l) + Beta(k) - sqrt(Sigma2(k));
                if ahigh - alow > 1 && Y(l) > alow && Y(l) < ahigh
                    Pxy(l, K + k) = 1 / (ahigh - alow);
                else
                    Pxy(l, K + k) = 0;
                end 
            else
                alow = Alpha(k + 1) * G(l) + Beta(k + 1);
                ahigh = Alpha(k) * G(l) + Beta(k);
                if Y(l) > alow && Y(l) < ahigh
                    Pxy(l, K + k) = exp( - (Y(l) - Alpha(k + K) * G(l) - Beta(k + K))^2 / 2 / Sigma2(k + K)) ...
                                / (sqrt(2 * pi * Sigma2(k + K)));
                else
                    Pxy(l, K + k) = 0;
                end
            end
        end
        for k = 1:K + K - 1
            Z1(l, k) = Pxy(l, k) * Pi(k) / sum(Pxy(l, :) .* Pi);
        end
    end

    % M-step: Update label marginal probabiliy and the regresssion parameters
    for k = 1:K * 2 - 1
        Pi(k) = sum(Z1(:, k)) / N;
        
        A = zeros(2);
        b = zeros(2, 1);
        A(1, 1) = sum(Z1(:, k) .* G2);
        A(1, 2) = sum(Z1(:, k) .* G);
        A(2, 1) = A(1, 2);
        A(2, 2) = sum(Z1(:, k));
        b(1) = sum(Z1(:, k) .* G .* Y);
        b(2) = sum(Z1(:, k) .* Y);

        if ~isnan(trace(A)) && ~isinf(trace(A)) && rank(A) == 2
        alpha_beta = A \ b;

        Alpha(k) = alpha_beta(1);
        Beta(k) = alpha_beta(2);
        if k > K && (Alpha(k) > Alpha(k-K) || Alpha(k) < Alpha(k-K+1))
        Alpha(k) = (Alpha(k-K)+Alpha(k-K+1))/2;
        Beta(k) = (Beta(k-K)+Beta(k-K+1))/2;
        end
        Sigma2(k) = sum( Z1(:, k) .* (Y - Alpha(k) * G - Beta(k)).^2 ) ...
                    / sum( Z1(:, k));
        end
    end
    Sigma2(Sigma2(1:K * 2 - 1) < 1) = 1;
    Sigma2(Sigma2(1:K * 2 - 1) > maxVar | isnan(Sigma2(1:K * 2 - 1))) = maxVar;
    
%     % DEBUG
%     R.Alpha = Alpha;
%     R.Beta = Beta;
%     R.Sigma2 = Sigma2;
%     R.Pi = Pi;
%     R.Q = Q(1:i);
%     R.Z1 = Z1;
%     R.X = X;
%     debug_segregress(X, Y, R);
%     pause(0.5);
end 

R.Alpha = Alpha(1:K);
R.Beta = Beta(1:K);
R.Sigma2 = Sigma2;
R.Pi = Pi;
R.Z1 = Z1;
R.AlphaC = Alpha;
R.BetaC = Beta;
R.SigmaC = sqrt(Sigma2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG
% debug_segregress(X, Y, R);
% figure(101), plot(R.Q);

