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

% Set average TOC and SW values for modeling
TOC_mean = mean(TOC);
SW_mean  = mean(SW);
Den_mean = mean(Den);

% Define porosity and crack density sample points
por = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121];
xp  = [0.15, 0.2, 0.28, 0.38, 0.5, 0.65, 0.8, 0.99];

%--------------------------------------------------------------------------
figure('Color','W')
xlabel('Shear Impedance');
ylabel('Vp/Vs Ratio');
zlabel('Density');
grid on;
view(45, 30);  

for i = 1:length(xp)
    for j = 1:length(por)
        [vp_pre(i,j), vs_pre(i,j), den(i,j)] = modeling( ...
            0.7, 0.05, 0.3, TOC_mean, SW_mean, por(j), Den_mean, xp(i));
    end
    Vw(i,:) = vp_pre(i,:) ./ vs_pre(i,:);        % Vp/Vs ratio
    Is(i,:) = vs_pre(i,:) .* den(i,:) / 1000;    % Shear impedance
    plot3(Is(i,:), Vw(i,:), den(i,:), 'r', 'LineWidth', 1);
    clearvars vp_pre vs_pre den;
    hold on
end
hold on

for i = 1:length(por)
    for j = 1:length(xp)
        [vp_pre(j,i), vs_pre(j,i), den(j,i)] = modeling( ...
            0.7, 0.05, 0.3, TOC_mean, SW_mean, por(i), Den_mean, xp(j));
    end
    Vw(:,i) = vp_pre(:,i) ./ vs_pre(:,i);        % Vp/Vs ratio
    Is(:,i) = vs_pre(:,i) .* den(:,i) / 1000;    % Shear impedance
    plot3(Is(:,i), Vw(:,i), den(:,i), 'k', 'LineWidth', 1);
    clearvars vp_pre vs_pre den;
    hold on
end
