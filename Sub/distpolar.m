function D2 = distpolar(ZI,ZJ)

global userPos0
ZI = ZI - userPos0;
ZJ = ZJ - userPos0;
m2 = length(ZJ);
D2 = zeros(m2, 1);
[thetai1, rhoi] = cart2pol(ZI(1), ZI(2));
thetai2 = atan2(ZI(3), rhoi);
rhoi = norm([rhoi ZI(3)]);
for i = 1:m2
    [thetaj1, rhoj] = cart2pol(ZJ(i, 1), ZJ(i, 2));
    thetaj2 = atan2(ZJ(i, 3), rhoj);
    theta1 = mod(thetai1 - thetaj1, 2 * pi);
    theta1 = min(theta1, 2 * pi - theta1) / pi * 180;
    theta2 = mod(thetai2 - thetaj2, 2 * pi);
    theta2 = min(theta2, 2 * pi - theta2) / pi * 180;
    rho = rhoi - norm([rhoj ZJ(i, 3)]);
    D2(i) = norm([theta1 theta2 rho]);
end
