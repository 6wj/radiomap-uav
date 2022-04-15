function L = likelihoodOf(k, h, Hs, Ps)

N = length(Ps);
N = min(N, size(Hs, 1));

L = 0;

for i = 1:N
    if ismember(k, Hs{i, 1})
        a = h(Hs{i, 1});
        b = Hs{i, 2};
        if sum(a(:) > b(:)) == 0
            L = L + Ps(i);
        else
            temp = h(k);
            h(k) = 0; % the k-th buiding height
            a = h(Hs{i, 1});
            if sum(a(:) > b(:)) == 0
                L = L + (1 - Ps(i));
            end
            h(k) = temp;
        end
    end
end