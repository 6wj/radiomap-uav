function [H, W] = optimizeH(R, Maps, H)
% Maps.BldPosMat and Maps.FolPosMat are abandoned here

K = size(H, 1);
M = size(H, 2);
N = size(R.X, 1);   % number of data samples
CofM = R.Hs; % Map cell number list of every measurement
MofC = R.S; % Measurement number list of every map cell
obj_all = R.Z1;

maxBldHeight = Maps.droneHeightMap;
weights = Maps.neighbourWeight;

kp = struct;
kp.kernelType = 'epanechnikov';
MAXLOOP = 10 * numel(H);      % maximum number of iterations

cnt = 1;
eps_h = 1;
H0 = H * 2;

while cnt < MAXLOOP && norm(H - H0, 'inf') > eps_h * 2
    H0 = H;
    for m = 1:M
        for k = K:-1:1
            hmax = maxBldHeight;
            hmin = 0;

            i_vec = MofC{m};
            if length(i_vec) < 2
                continue;
            end
            
            Vk_vec = zeros(1, 5);
            H(k, m) = hmin;
            Vk_vec(1) = Val(min(maxBldHeight, max(0, H)), ...
                CofM(i_vec, :), obj_all(i_vec, :), weights);
            H(k, m) = hmax;
            Vk_vec(5) = Val(min(maxBldHeight, max(0, H)), ...
                CofM(i_vec, :), obj_all(i_vec, :), weights);
            while hmax - hmin > eps_h
                hs = hmin: (hmax - hmin) / 4:hmax;
                for ih = 2:4
                    H(k, m) = hs(ih);
                    Vk_vec(ih) = Val(min(maxBldHeight, max(0, H)), ...
                        CofM(i_vec, :), obj_all(i_vec, :), weights);
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
            hmk = hs(Imax);
            H(k, m) = hmk;

            cnt = cnt + 1;
        end
    end
end

% Generate w, LoS/NLoS labels from heights
neighbourSize = size(CofM, 2) / 2;
W = zeros(N, 2*(K+1));
for i = 1:N
    for k = 0:K
        for nr = 1:neighbourSize
            if sum(vec(H(k+1:K, CofM{i, nr * 2 - 1}) > ...
                    repmat(CofM{i, nr * 2}, K-k, 1))) == 0 && ...
                    (k == 0 || ...
                    sum(H(k, CofM{i, nr * 2 - 1}) > CofM{i, nr * 2}) > 0)
                W(i, 2*k+1) = W(i, 2*k+1) + weights(nr);
            end
        end
    end
end
W(:, 2:2:end) = W(:, 1:2:end);

end

function V = Val(h, Hs, Ps, weights)

N = size(Ps, 1);
N = min(N, size(Hs, 1));
K = size(h, 1);
neighbourSize = size(Hs, 2) / 2;

V = 0;

for i = 1:N
    for k = 0:K
        for nr = 1:neighbourSize
            if sum(vec(h(k+1:K, Hs{i, nr * 2 - 1}) > ...
                    repmat(Hs{i, nr * 2}, K-k, 1))) == 0 && ...
                    (k == 0 || ...
                    sum(h(k, Hs{i, nr * 2 - 1}) > Hs{i, nr * 2}) > 0)
                V = V + weights(nr) * Ps(i, k+1);
            end
        end
    end
end
if N == 0, V = - intmax; end

end
