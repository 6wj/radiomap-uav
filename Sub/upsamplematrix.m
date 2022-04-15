function B = upsamplematrix(A, dim)

N = dim(1);
M = dim(2);

[n, m] = size(A);
B = zeros(N, M);
for i = 1:N
    for j = 1:M
        B(i, j) = A( max(1, round(i * n / N)), ...
                     max(1, round(j * m / M)) );
    end
end