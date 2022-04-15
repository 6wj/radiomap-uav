clear, close all

DATA = load('radiomap_shanghai100tx_2.5GHz.mat');
% DATA = load('radiomap_shanghai115tx_28GHz.mat');

PosUE = DATA.PosUE;
droneHeightSample = 50;
droneHeightMap = 50;
nUser = 30;
nDrone = 30 * length(droneHeightSample);
nMeas = nUser * nDrone;
nDrone = round(nDrone / length(droneHeightSample));
noiseGuass = 3; % dB std
Xdata = zeros(nMeas, 6); 
Gdata = zeros(nMeas, 1);
Ddata = zeros(nMeas, 1);
Gtrue = cell(nUser, 1);
cnt = 0;

figure(1);
userPosAll = zeros(size(PosUE));
userids = sort(randperm(length(PosUE), nUser));
for id = 1:length(userids)
    % i must increase
    i = userids(id);
    % Training data
    for j = 1:length(droneHeightSample)
        DroneHeight = droneHeightSample(j);
        x = PosUE(i,1); y = PosUE(i,2); z = PosUE(i,3);
        I = (DATA.RadioMap(:, 6) == DroneHeight) ...
            & (DATA.RadioMap(:, 1) == x) ...
            & (DATA.RadioMap(:, 2) == y)...
            & (DATA.RadioMap(:, 3) == z);
        Rm2D = DATA.RadioMap(I, :);
        Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
        [Xmat,Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
        Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
        UserPos = [x-min(Xvec)+1 y-min(Yvec)+1 z];
        Xmat = Xmat - min(Xvec) + 1;
        Ymat = Ymat - min(Yvec) + 1;
        % Data pre-processing
        Zmat(isnan(Zmat)) = min(min(Zmat));
        Zmat = Zmat + randn(size(Zmat)) * noiseGuass;
        Xvec = Xmat(:);
        Yvec = Ymat(:);
        Zvec = Zmat(:);
        dist = log10(vecnorm([Xvec Yvec DroneHeight*ones(length(Xvec),1)] - ...
            ones(length(Xvec),3).*UserPos, 2, 2));
        UserPos = [ceil(UserPos(1:2)) z];
        if ~exist('idx', 'var')
            idx = randperm(length(Xvec), nDrone);
        end
        DronePositon = ceil([Xvec(idx) Yvec(idx)]);
        DroneNum = length(DronePositon);
        cnt = cnt + DroneNum;
        Xdata(cnt-DroneNum+1:cnt, :) = [DronePositon ...
            DroneHeight.*ones(DroneNum,1) ones(DroneNum,3).*UserPos];
        Gdata(cnt-DroneNum+1:cnt, :) = Zvec(idx);
        Ddata(cnt-DroneNum+1:cnt, :) = dist(idx);
        userPosAll(i, :) = UserPos;
    end
    % Test data
    DroneHeight = droneHeightMap;
    x = PosUE(i,1); y = PosUE(i,2); z = PosUE(i,3);
    I = (DATA.RadioMap(:, 6) == DroneHeight) ...
        & (DATA.RadioMap(:, 1) == x) ...
        & (DATA.RadioMap(:, 2) == y)...
        & (DATA.RadioMap(:, 3) == z);
    Rm2D = DATA.RadioMap(I, :);
    Xvec = Rm2D(:, 4); Yvec = Rm2D(:, 5); Zvec = Rm2D(:, end);
    [Xmat, Ymat] = meshgrid((min(Xvec):5:max(Xvec)),(min(Yvec):5:max(Yvec)));
    Zmat = griddata(Xvec, Yvec, Zvec, Xmat, Ymat);
    UserPos = [x-min(Xvec)+1 y-min(Yvec)+1 z];
    Xmat = Xmat - min(Xvec) + 1;
    Ymat = Ymat - min(Yvec) + 1;
    % Data pre-processing
    Xvec = Xmat(:);
    Yvec = Ymat(:);
    Zvec = Zmat(:);
    Zmat(isnan(Zmat)) = min(Zvec);
    Gtrue{i} = Zmat;
    figure(1);
    surf(Xmat, Ymat, Zmat, 'edgecolor', 'none', 'FaceColor', 'interp');
    view(0, 90);
    axis square
    xlim([0, max(Xvec)]);
    ylim([0, max(Yvec)]);
    pause(1);
end
