function [covBlds, covBldZs] = covBldZ(DronePos, UePos, lenX, noZ)
% 
% OUTPUT
%   covBlds         Index set of buildings under the line segment
%   covBldZs        Z values on the 3D line segment 

xu = UePos(1);
yu = UePos(2);
zu = UePos(3);

xd = DronePos(1);
yd = DronePos(2);
zd = DronePos(3);

% find the buildings under the line
covBlds0 = zeros(1, abs(xu - xd) + abs(yu - yd));
cnt = 0;
if xd < xu
    I = xd:xu;
elseif xu < xd
    I = xu:xd;
else
    I = [];
end
for xc = I
    yc = round((yd - yu) * (xc - xu) / (xd - xu) + yu);
    cnt = cnt + 1;
    covBlds0(cnt) = (yc - 1) * lenX + xc;
end

if yd < yu
    I = yd:yu;
elseif yu < yd
    I = yu:yd;
else
    I = [];
end
for yc = I
    xc = round((yc - yu) * (xd - xu) / (yd - yu) + xu);
    cnt = cnt + 1;
    covBlds0(cnt) = (yc - 1) * lenX + xc;
end
covBlds = union(covBlds0, []);  % Indices of all the buildings under the line

if nargin > 3 && noZ, return, end
ncovBlds = length(covBlds);
covBldZs = zeros(1, ncovBlds);

for i = 1:ncovBlds
    j = covBlds(i);
    
    yb = floor((j - 1) / lenX) + 1;
    xb = j - (yb - 1) * lenX;
    zb = zu + (zd - zu) * norm([xu yu] - [xb yb], 2) / norm([xu yu] - [xd yd], 2);
    
    covBldZs(i) = zb;
end

                  