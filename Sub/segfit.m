function R = segfit(X, Y, debug)
% Segmented linearly fit the data to two linear models
%
% Joint data clustering and regression based on EM algorithms with
% parameter controls (from heuristics of the radio map application). 
% 
% We cluster the data to the model: y = a_i * x + b_i, where a_i and b_i
% are model parameters for i = {1,2}. The intuition obtained from analyzing
% real data suggest motivates the design of this script. The method is only
% robust for just two clusters. 
% 
% Input: X and Y are vectors.
% Output: R.A = [a1, a2], R.B = [b1, b2] are the model parameters; R.C = is
% 	a vector of soft label where c = 1 means cluster 1, c = 0 means cluster
% 	2; R.S = [s1, s2] is the standard deviation for the two clusters.
%
% Modified from segregressRM.m 
% Date: March 20, 2020.

maxSigma1 = inf;   % The STD of cluster 1 should not be lower than this value
DEBUG_MODE = 0;

if nargin < 1
    DEBUG_MODE = 1;
    N = 10000;                     % Only analyze a subset of data
    DATA = load('radiomap.mat'); 
    nRow = size(DATA.RadioMap, 1); N = min(nRow, N);
    I = randperm(nRow); I = I(1:N);
    Rm = DATA.RadioMap(I, :);
    X = log10(Rm(:, 7)); Y = Rm(:, end); K = 2;
elseif nargin > 2
    DEBUG_MODE = debug;
end

N = length(X);

% Initialization from a global fit
P = polyfit(X, Y, 1);
Z = Y - (P(1) * X + P(2));
I1 = Z > 0; I2 = Z < 0;
% Sigma2_1 = var(Z(I1)); Sigma2_2 = var(Z(I2));

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

F = F(1:cnt);

if DEBUG_MODE
    Ap = [P1(1), P2(1)]; Bp = [P1(2), P2(2)];

    figure(4000), plot(F, 'linewidth', 2); 
    title('Convergence'); xlabel('Iteration'); ylabel('Objective (distance)');

    figure(4001), plot(X(I1), Y(I1), '.', 'markersize', 1);hold on
    plot(X(I2), Y(I2), '.', 'markersize', 1);
    t = min(X):(max(X) - min(X))/10:max(X); cmap = colormap('prism');
    for k = 1:2
        plot(t, Ap(k) * t + Bp(k), 'color', cmap(k, :));
    end
    hold off
    xlabel('log-distance'); ylabel('path gain');
    title(sprintf('a = (%3.2f, %3.2f), b = (%3.2f, %3.2f), s = (%3.2f, %3.2f)', ...
        Ap(1), Ap(2), Bp(1), Bp(2), sqrt(S1), sqrt(S2)));
end

Ap = [P1(1), P2(1)]; Bp = [P1(2), P2(2)];

R.Alpha = Ap;
R.Beta = Bp;
R.Sigma = sqrt([S1, S2]);
R.Z1 = C1;
