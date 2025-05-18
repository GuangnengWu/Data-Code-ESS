clear; clc; tic;

load('LogDataALL.mat')

% Load logging data
VP    = LogDataALL(1:900, 1);    % P-wave velocity (m/s)
VS    = LogDataALL(1:900, 2);    % S-wave velocity (m/s)
Den   = LogDataALL(1:900, 3);    % Density (g/cm³)
Vli   = LogDataALL(1:900, 4);    % Calcite volume fraction
POR   = LogDataALL(1:900, 5);    % Porosity
XP    = LogDataALL(1:900, 6);    % Pore aspect ratio

% Set average TOC and SW values for modeling
TOC_mean = 0.025;     % Average Total Organic Carbon
SW_mean  = 0.96;      % Average water saturation
Den_mean = mean(Den); % Average density

% Calculate elastic attributes
VRobs = VP ./ VS;          % Vp/Vs ratio
ISob  = VS .* Den ./ 1000; % Shear impedance (scaled)

% Define parameter ranges
vli=0.005:0.005:0.70;
por=0.001:0.005:0.121;
xp=0.05:0.01:1;

% Run forward modeling to generate rock physics templates
for k = 1:length(vli)
    for j = 1:length(por)
        for i = 1:length(xp)
            vsh = 1 - vli(k) - 0.075;
            vqu = 1 - vsh - vli(k) - 0.025;
            lime(i,j,k)     = vli(k);
            porosity(i,j,k) = por(j);
            poar(i,j,k)     = xp(i);
            [vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(j), Den_mean, xp(i));
        end
        Vw(:,j,k) = vp(:,j,k) ./ vs(:,j,k); 
        Is(:,j,k) = vs(:,j,k) .* rho(:,j,k) ./ 1000;   
        hold on
    end
end

% Flatten the synthetic elastic parameters to 1D vectors
x = Is;   xx = x(:);
y = Vw;   yy = y(:);
z = rho;  zz = z(:);

% Flatten the modeled physical properties
aa = lime;     a = aa(:);
bb = porosity; b = bb(:);
cc = poar;     c = cc(:);

% Create interpolants to map elastic attributes to physical properties
F1 = scatteredInterpolant(xx, yy, zz ,a, 'linear', 'linear');
F2 = scatteredInterpolant(xx, yy, zz ,b, 'linear', 'linear');
F3 = scatteredInterpolant(xx, yy, zz ,c, 'linear', 'linear');

% Predict physical properties using the interpolants
vlipre = F1(ISob , VRobs , Den);vlipre(vlipre < 0.001) = 0.001;
porpre = F2(ISob , VRobs , Den);porpre(porpre < 0.001) = 0.001;
xppre  = F3(ISob , VRobs , Den);xppre(xppre < 0.001) = 0.001;  xppre(xppre > 1) = 1;

% Visualization
figure('Color','W')
depth = 1:900;

subplot(1,5,2)
plot(Vli, depth, 'r', vlipre, depth, 'k'); 
set(gca, 'YDir', 'reverse');
xlim([0 1]);
xlabel('Calcite Volume Fraction (Vli)');
ylabel('Sample Index');
legend('Measured', 'Predicted');
title('Calcite Volume Fraction');

subplot(1,5,3)
plot(POR, depth, 'r', porpre, depth, 'k'); 
set(gca, 'YDir', 'reverse');
xlabel('Porosity');
legend('Measured', 'Predicted');
title('Porosity');

subplot(1,5,4)
plot(XP, depth, 'r', xppre, depth, 'k'); 
set(gca, 'YDir', 'reverse');
xlabel('Aspect Ratio');
legend('Measured', 'Predicted');
title('Aspect Ratio');