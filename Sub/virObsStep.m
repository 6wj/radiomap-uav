function [BldMapHat, w] = virObsStep(R, Maps, H0)
% Building/Foliage height estimation 
% [BldMapHat, FolMapHat] = virtualBldMap(R, Maps, Map0)
% 
% INPUT
%   R           Output from the EM algorithm 
%   Maps        A structure that contains: 
%       Maps.BldPosMat       An indicator matrix for the positions of buildings
%       Maps.FolPosMat       An indicator matrix for the positions of trees
%       Maps.meterPerPixel   A scalar for how many meters for a matrix
%                            element
%   H0          Initial building height values
 
N = size(R.X, 1);   % number of data samples
Hs = R.Hs;
S = R.S;

% - Building height
if exist('H0', 'var')
    maxBldHeight = max(50, max(vec(H0)));  % meter
else
    maxBldHeight = 50;
end
BldPos = find(Maps.BldPosMat > 0);
nBldPos = length(BldPos);
if nargin < 3
    B0 = ones(1, nBldPos) * maxBldHeight;
else 
    B0 = H0(Maps.BldPosMat > 0);
    B0(B0 > maxBldHeight) = maxBldHeight;
end

%
Hhat = bldHeight3(Hs, R.Z1, S, maxBldHeight, B0, 0);
BldMapHat = zeros(size(Maps.BldPosMat));
BldMapHat(BldPos) = Hhat;

% Generate w, LoS/NLoS labels from heights
neighbourSize = size(Hs, 2) / 2;
w = zeros(N, 1);
for i = 1:N
    losInd = 0;
    for nr = 1:neighbourSize
        losInd = losInd + (isempty(Hs{i, nr * 2 - 1}) | ...
            sum(Hhat(Hs{i, nr * 2 - 1}) > Hs{i, nr * 2}') == 0);
    end
    w(i) = losInd / neighbourSize;
end

end

function Hhat = bldHeight3(Hs, Ps, S, maxBldHeight, H0, show_figures)
% Version 3: Added initial value H0
% Version 2: accelerate the algorithm by introducing S (the set of related 
% measurements for each building). Moreover, Hs is now a cell array, for 
% row, it contains the indices of the buildings and the Z values
%
% Building height estimation
%
% INPUT     
%   Hs              A matrix (# of data samples) * (# of buildings),
%                   where the i-th row gives the height vector h_i. Recall
%                   the joint CDF is expressed as 
%                           P_i = Prob{h <= h_i}
%
%   Ps              A vector with (# of data samples) elements,
%                   each i-th entry representing the probability for h_i 
%
%   maxBldHeight    Maximum building height 
% 
%   show_figures    (Optional) {0, 1}
%   

if nargin < 4
    show_figures = 0;
end

if iscell(Hs)
    N = size(Ps, 1);
    Nbld = length(H0(:));
else
    [N, Nbld] = size(Hs);
end

kp = struct;
kp.kernelType = 'epanechnikov';
eps_h = 1;       % Buildin height (ideal) precision in meters
MAXLOOP = 10 * Nbld;      % maximum number of iterations
    
% H_array = zeros(MAXLOOP, Nbld);

% H0 = maxBldHeight * ones(1, Nbld);
H = H0;

cnt = 1;
% H_array(cnt, :) = H;
H0 = H * 2;

while cnt < MAXLOOP && norm(H0 - H, 'inf') > eps_h * 1.5
    if show_figures
        fprintf('Building height estimation iteration: %2d\nBld 0000', ceil((cnt + 1) / Nbld));
    end
    
    H0 = H;
    for k = 1:Nbld
        if show_figures
            fprintf('\b\b\b\b%4d', k);
        end
        Hk = H;
        hmax = maxBldHeight;
        hmin = 0;

        kI = S{k};
        if length(kI) < 2
            continue;
        end
        
        Vk_vec = zeros(1, 5);   % Only focus on the value that invovles the kth building
        Hk(k) = hmin;
        Vk_vec(1) = Val(min(maxBldHeight, max(0, Hk)), Hs(kI, :), Ps(kI, :));
        Hk(k) = hmax;
        Vk_vec(5) = Val(min(maxBldHeight, max(0, Hk)), Hs(kI, :), Ps(kI, :));
        while hmax - hmin > eps_h
            hs = hmin: (hmax - hmin) / 4:hmax;
            for ih = 2:4
                Hk(k) = hs(ih);
                Vk_vec(ih) = Val(min(maxBldHeight, max(0, Hk)), Hs(kI, :), Ps(kI, :));
%{
                Xtr = [min(maxBldHeight, max(0, Hk(k) - 0.9));
                    min(maxBldHeight, max(0, Hk(k) - 0.5));
                    min(maxBldHeight, max(0, Hk(k)));
                    min(maxBldHeight, max(0, Hk(k) + 0.5));
                    min(maxBldHeight, max(0, Hk(k) + 0.9))];
                Ytr = [
                    Val(min(maxBldHeight, max(0, Hk - 0.9)), Hs(kDataIds, :), Ps(kDataIds, :));
                    Val(min(maxBldHeight, max(0, Hk - 0.5)), Hs(kDataIds, :), Ps(kDataIds, :));
                    Val(min(maxBldHeight, max(0, Hk)), Hs(kDataIds, :), Ps(kDataIds, :));
                    Val(min(maxBldHeight, max(0, Hk + 0.5)), Hs(kDataIds, :), Ps(kDataIds, :));
                    Val(min(maxBldHeight, max(0, Hk + 0.9)), Hs(kDataIds, :), Ps(kDataIds, :))];
                Vk_vec(ih) = localPolyRegression(min(maxBldHeight, max(0, Hk(k))),Xtr,Ytr,'',1,kp);
%}
            end
            [maxVh, ~] = max(Vk_vec);
            % Option 1: tends to over estiamte the building height
%             Imax = find(abs(Vk_vec - maxVh) < 1e-9, 1, 'last');
            % Option 2: guess the medium value
            Imax0 = find(abs(Vk_vec - maxVh) < 1e-9);
            Imax = ceil(mean(Imax0));
            if Imax == 3
                i_left = 2;
                i_right = 4;
            elseif Imax < 3
                i_left = 1;
                i_right = 3;
            else
                i_left = 3;
                i_right = 5;
            end
               
            hmin = hs(i_left);
            hmax = hs(i_right);
            
            Vk_vec0 = Vk_vec;
            Vk_vec(1) = Vk_vec0(i_left);
            Vk_vec(5) = Vk_vec0(i_right);
        end
        hk = hs(Imax);
        Hk(k) = hk;
        H = Hk;
        
        cnt = cnt + 1;
%         H_array(cnt, :) = H;
       
    end
    if show_figures
        fprintf('\n');
    end

%     H_array(cnt, :) = H;
    
end
% Hhat = H_array(cnt, :);
Hhat = H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization for 2 building case ---------------------------------------
% Illustration: Empirical CDF for two building case
if show_figures
    if Nbld == 2
        figure(5),
        cmap = jet(100);
        for i = 1:N
            if Ps(i) > 0.5
                plot3(Hs(i, 1), Hs(i, 2), Ps(i), 'o', 'color', ...
                            cmap(floor((max(0, Ps(i) - 1e-9)) * 100) + 1, :));
            else
                plot3(Hs(i, 1), Hs(i, 2), Ps(i), 'x', 'color', ...
                            cmap(floor((max(0, Ps(i) - 1e-9)) * 100) + 1, :));
            end
            hold on
        end
        xlim([0, maxBldHeight]);
        ylim([0, maxBldHeight]);
        set(gca, 'FontSize', 14);
        set(gca, 'XTick', 0:10:50);
        set(gca, 'YTick', 0:10:50);
        xlabel('h1');
        ylabel('h2');
        view([0 90]);
        axis square
        hold off
        
        hold on
        plot3(H_array(1:cnt, 1), H_array(1:cnt, 2), ones(cnt, 1), 'gs-');
        hold off
    
    end

    % Illustration: Objective function for two building cass
    if Nbld == 2
        figure(6),
        st = 50;
        [H1, H2] = meshgrid( (1/st:1/st:1) * maxBldHeight, (1/st:1/st:1) * maxBldHeight);
        nH = size(H1, 1);
        Vmat = zeros(nH);
        for i = 1:nH
            for j = 1:nH
                h = [H1(i, j), H2(i, j)];
                Vmat(i, j) = Val(h, Hs, Ps);
            end
        end
        surf(H1, H2, Vmat, 'EdgeColor', 'none');
        xlabel('h1 (m)');
        ylabel('h2 (m)');
        % zlabel('Objective value');
        title('Objective value');
        axis square
        view(0, 90);
    	
        hold on
        plot3(H_array(1:cnt, 1), H_array(1:cnt, 2), max(Vmat(:)) * ones(cnt, 1), 'ks-');
        hold off
    
    end
end

end

function V = Val(h, Hs, Ps)

N = size(Ps, 1);
N = min(N, size(Hs, 1));
neighbourSize = size(Hs, 2) / 2;

V = 0;

for i = 1:N
    losInd = 0;
    for nr = 1:neighbourSize
        if sum(h(Hs{i, nr * 2 - 1}) > Hs{i, nr * 2}') == 0
            losInd = losInd + 1;
        end
    end
    losInd = losInd / neighbourSize;
    V = V + losInd * Ps(i, 1) + (1 - losInd) * Ps(i, 2);
end
if N == 0, V = - intmax; end

end
