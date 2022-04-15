function V = Val(h, Hs, Ps)
%

N = size(Ps, 1);
N = min(N, size(Hs, 1));
M = size(Ps, 2);

V = 0;

if M == 0, warning('M = 0'), end
if iscell(Hs) && M == 2
    for i = 1:N
        a = h(Hs{i, 1});
        b = Hs{i, 2};
        if sum(a(:) > b(:)) == 0
            V = V + Ps(i, 1);
        else
            V = V + Ps(i, 2);
        end
    end
    if N == 0, V = - intmax; end
    
elseif iscell(Hs)
    for i = 1:N
        a = h(Hs{i, 1});
        b = Hs{i, 2};
        if sum(a(:) > b(:)) == 0
            V = V + Ps(i);
        else
            V = V + (1 - Ps(i));
        end
    end
    
else
    for i = 1:N
        if sum(h(Hs(i, :) > 0) > Hs(i, Hs(i, :) > 0)) == 0
            V = V + Ps(i);
        else
            V = V + (1 - Ps(i));
        end
    end
    
end
