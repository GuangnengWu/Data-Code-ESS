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

por=[0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.121];
vli=[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.71];

%--------------------------------------------------------------------------
figure('Color','W')
xlabel('Shear Impedance');
ylabel('Vp/Vs Ratio');
zlabel('Density');
grid on;
view(45, 30);  
for i = 1:length(vli)
    for j = 1:length(por)
        vsh = 1 - vli(i) - 0.075;
        vqu = 1 - vsh - vli(i) - 0.025;
        [vp_pre(i,j), vs_pre(i,j), den_mat(i,j)] = modeling( ...
            vsh, vqu, vli(i), TOC_mean, SW_mean, por(j), Den_mean, 0.25);
    end
    Vw = vp_pre(i,:) ./ vs_pre(i,:);        % Vp/Vs ratio
    Is = vs_pre(i,:) .* den_mat(i,:) / 1000; % Shear impedance
    plot3(Is, Vw, den_mat(i,:), 'r', 'LineWidth', 1);
    hold on
end

for i = 1:length(por)
    for j = 1:length(vli)
        vsh = 1 - vli(j) - 0.075;
        vqu = 1 - vsh - vli(j) - 0.025;
        [vp_pre(j,i), vs_pre(j,i), den_mat(j,i)] = modeling( ...
            vsh, vqu, vli(j), TOC_mean, SW_mean, por(i), Den_mean, 0.25);
    end
    Vw = vp_pre(:,i) ./ vs_pre(:,i);        % Vp/Vs ratio
    Is = vs_pre(:,i) .* den_mat(:,i) / 1000; % Shear impedance
    plot3(Is, Vw, den_mat(:,i), 'k', 'LineWidth', 1);
    hold on
end