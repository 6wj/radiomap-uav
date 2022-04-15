function H = mlHeight(k, h, Hs, Ps, S, ds, maxHeight, nBldBlob, blobsLabel)

%global metrics;

if nargin == 9
    % Building lines belong to the same building blob
    bldBlob = nBldBlob{blobsLabel(k), 1};
    h_compare = h;
    h_compare(bldBlob) = 0;
else
    h_compare = h;
    h_compare(k) = 0;
end

% Generate feasible set of the k-th height
L = containers.Map('KeyType', 'double', 'ValueType', 'double');
for j = 1:length(S{k})
    i = S{k}(j);
    if mod(i, ds) ~= 0
        continue;
    end
    idOfK = find(k == Hs{i, 1});
    L(Hs{i, 2}(idOfK(1))) = 0;
end
heights_vector = cell2mat(keys(L)); % All possible heights

cnt = 0;
hmin = 1;
hmax = length(heights_vector);
L_vector = zeros(1, 5);
hs = floor(hmin:(hmax-hmin)/4:hmax);
hs_loop = hs;

while hmax - hmin > 3
    for ih = 1:length(hs_loop)
        cnt = cnt + 1;
        for j = 1:length(S{k})
            i = S{k}(j);
            if mod(i, ds) ~= 0
                continue;
            end
            idOfK = find(k == Hs{i, 1});

            if sum(h_compare(Hs{i, 1}) > Hs{i, 2}) == 0
                %metrics(i) = 1;
                if heights_vector(hs_loop(ih)) <= Hs{i, 2}(idOfK(1)) % L_LOS
                    L_vector(cnt) = L_vector(cnt) + ...
                        Ps(i, 1); %((Ps(i, 1) + Ps(i, 2) * (nargin == 9)) > 0.5);
                else % L_NLOS
                    L_vector(cnt) = L_vector(cnt) + ...
                        1 - Ps(i, 1); %((1 - Ps(i, 1) - Ps(i, 2) * (nargin == 9)) > 0.5);
                end
            end
        end
    end
    [~, idOfMax] = max(L_vector);
    if idOfMax == 3
        i_left = 2;
        i_right = 4;
    elseif idOfMax < 3
        i_left = 1;
        i_right = 3;
    else
        i_left = 3;
        i_right = 5;
    end
    hmin = hs(i_left);
    hmax = hs(i_right);
    hs = floor(hmin:(hmax-hmin)/4:hmax);
    hs_loop = hs(2:4);
    L_vector_temp = L_vector;
    L_vector(1) = L_vector_temp(i_left);
    L_vector(5) = L_vector_temp(i_right);
    L_vector(2:4) = 0;
    cnt = 1;
end

if (nargin ~= 9 && length(hs) >= 3 && hs(3) > 0) || ... % For foliage height
        hmax >= 3 && heights_vector(hs(3)) > 0 && heights_vector(hs(3)) < maxHeight
    H = min(maxHeight, max(0, heights_vector(hs(3))));
elseif nargin == 9
    H = mean(h(bldBlob));
    H = H * (H <= maxHeight);
else
    H = h(k);
end