clear; clc; tic;

load('LogData1.mat')  % Load logging data

% Extract the first 561 samples of well log data
Depth = LogData(1:561, 1);    % Depth (m)
VP    = LogData(1:561, 2);    % P-wave velocity (m/s)
VS    = LogData(1:561, 3);    % S-wave velocity (m/s)
Den   = LogData(1:561, 4);    % Density (g/cm³)
Vli   = LogData(1:561, 5);    % Calcite content
POR   = LogData(1:561, 6);    % Porosity
Vqu   = LogData(1:561, 7);    % Quartz content
Vsh   = LogData(1:561, 8);    % Clay content
SW    = LogData(1:561, 9);    % Water saturation
TOC   = LogData(1:561,10);    % Total Organic Carbon (TOC) content

% Set average values
TOC_mean = mean(TOC);
SW_mean  = mean(SW);
Den_mean = mean(Den);  % Use mean density in modeling

% Define parameter grid
xp  = [0.08, 0.095, 0.12, 0.15, 0.2, 0.35, 0.99];  % Fracture
vli = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.71];  % Calcite

% Initialize figure
figure('Color','W')
xlabel('Shear Impedance');
ylabel('Vp/Vs Ratio');
zlabel('Density');
grid on;
view(45, 30);  

for i = 1:length(vli)
    for j = 1:length(xp)
        vsh = 1 - vli(i) - 0.075;
        vqu = 1 - vsh - vli(i) - 0.025;
        [vp_mat(i,j), vs_mat(i,j), den_mat(i,j)] = modeling( ...
            vsh, vqu, vli(i), TOC_mean, SW_mean, 0.03, Den_mean, xp(j));
    end
    VpVs = vp_mat(i,:) ./ vs_mat(i,:);
    Is   = vs_mat(i,:) .* den_mat(i,:) / 1000;
    plot3(Is, VpVs, den_mat(i,:), 'r', 'LineWidth', 1);
    hold on
end

for i = 1:length(xp)
    for j = 1:length(vli)
        vsh = 1 - vli(j) - 0.075;
        vqu = 1 - vsh - vli(j) - 0.025;
        [vp_mat(j,i), vs_mat(j,i), den_mat(j,i)] = modeling( ...
            vsh, vqu, vli(j), TOC_mean, SW_mean, 0.03, Den_mean, xp(i));
    end
    VpVs = vp_mat(:,i) ./ vs_mat(:,i);
    Is   = vs_mat(:,i) .* den_mat(:,i) / 1000;
    plot3(Is, VpVs, den_mat(:,i), 'k', 'LineWidth', 1);
    hold on
end