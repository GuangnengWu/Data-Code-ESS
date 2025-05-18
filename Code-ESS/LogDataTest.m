clear; clc; tic;

load('LogData1.mat')  % Load logging data

% Extract well log data (first 561 points)
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

% Plot logging curves
figure
subplot(1,9,1)
plot(VP, Depth, 'r', 'LineWidth', 1)
xlim([2500 4500]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('VP');

subplot(1,9,2)
plot(VS, Depth, 'r', 'LineWidth', 1)
xlim([1500 2500]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('VS');

subplot(1,9,3)
plot(Den, Depth, 'r', 'LineWidth', 1);
xlim([2.4 2.7]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('Density');

subplot(1,9,4)
plot(Vsh, Depth, 'b', 'LineWidth', 1);
xlim([0 1]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('Clay');

subplot(1,9,5)
plot(Vli, Depth, 'b', 'LineWidth', 1);
xlim([0 1]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('Calcite');

subplot(1,9,6)
plot(smooth(Vqu), Depth, 'b', 'LineWidth', 1);
xlim([0 0.2]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('Quartz');

subplot(1,9,7)
plot(TOC, Depth, 'b', 'LineWidth', 1);
xlim([0 0.05]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('TOC');

subplot(1,9,8)
plot(POR, Depth, 'b', 'LineWidth', 1);
xlim([0 0.15]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('Porosity');

subplot(1,9,9)
plot(SW, Depth, 'b', 'LineWidth', 1);
xlim([0 1]); ylim([3530 3600]); set(gca, 'YDir', 'reverse');
xlabel('SW');









