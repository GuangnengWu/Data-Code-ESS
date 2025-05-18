clear; clc; tic;

load('LogDataALL.mat')

VP    = LogDataALL(1:900, 1);    % P-wave velocity (m/s)
VS    = LogDataALL(1:900, 2);    % S-wave velocity (m/s)
Den   = LogDataALL(1:900, 3);    % Density (g/cm³)
Vli   = LogDataALL(1:900, 4);    % Calcite volume fraction
POR   = LogDataALL(1:900, 5);    % Porosity
XP    = LogDataALL(1:900, 6);    % Pore aspect ratio

% Set average TOC and SW values for modeling
TOC_mean = 0.025;    % Average TOC
SW_mean  = 0.96;     % Average water saturation
Den_mean = mean(Den); % Average density

% Calculate elastic attributes
VRobs = VP ./ VS;          % Vp/Vs ratio
ISob  = VS .* Den ./ 1000; % Shear impedance (converted to proper scale)

% Define parameter ranges
por = [0.001, 0.04, 0.08, 0.12, 0.15];
xp  = [0.18, 0.28, 0.4, 0.6, 1];
vli = [0.01, 0.25, 0.45, 0.6, 0.7];

%% Plot modeled rock physics template lines
figure('Color','w'); 
grid on;
scatter3(ISob, VRobs, Den, 20, POR, ...
    'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat')
xlabel('Shear Impedance'); ylabel('Vp/Vs Ratio'); zlabel('Density');
xlim([3 7]); ylim([1.2 2]); zlim([2.3 2.7]);
view(-40, 24);colormap(jet);colorbar;hold on        
for k=1:length(vli)
    for j=1:length(por)
        for i=1:length(xp)
vsh=1-vli(k)-0.075;
vqu=1-vsh-vli(k)-0.025;
[vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(j), Den_mean, xp(i));
        end
        Vw(:,j,k)=vp(:,j,k)./vs(:,j,k); 
        Is(:,j,k)=vs(:,j,k).*rho(:,j,k)./1000;   
        plot3(Is(:,j,k),Vw(:,j,k),rho(:,j,k), 'k','LineWidth', 1);
        hold on
    end
hold on
    for j=1:length(xp)
        for i=1:length(por)
vsh=1-vli(k)-0.075;
vqu=1-vsh-vli(k)-0.025;
vqu(vqu<0)=0.001;
[vp(j,i,k), vs(j,i,k), rho(j,i,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(i), Den_mean, xp(j));
        end
        Vw(j,:,k)=vp(j,:,k)./vs(j,:,k); 
        Is(j,:,k)=vs(j,:,k).*rho(j,:,k)./1000;   
        plot3(Is(j,:,k),Vw(j,:,k),rho(j,:,k), 'k', 'LineWidth',1);
        hold on
    end
end
hold on
for k=1:length(por)
    for j=1:length(vli)
        for i=1:length(xp)
vsh=1-vli(j)-0.075;
vqu=1-vsh-vli(j)-0.025;
vqu(vqu<0)=0.001;
[vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(j), TOC_mean, SW_mean, por(k), Den_mean, xp(i));
        end
        Vw(:,j,k)=vp(:,j,k)./vs(:,j,k); 
        Is(:,j,k)=vs(:,j,k).*rho(:,j,k)./1000;    
        plot3(Is(:,j,k),Vw(:,j,k),rho(:,j,k), 'k', 'LineWidth', 1);
        clearvars vp vs rho Vw Is;
        hold on
    end
hold on
    for j=1:length(xp)
        for i=1:length(vli)
vsh=1-vli(i)-0.075;
vqu=1-vsh-vli(i)-0.025;
vqu(vqu<0)=0.001;
[vp(j,i,k), vs(j,i,k), rho(j,i,k)] = modeling(vsh, vqu, vli(i), TOC_mean, SW_mean, por(k), Den_mean, xp(j));
        end
        Vw(j,:,k)=vp(j,:,k)./vs(j,:,k);
        Is(j,:,k)=vs(j,:,k).*rho(j,:,k)./1000;  
        plot3(Is(j,:,k),Vw(j,:,k),rho(j,:,k), 'k', 'LineWidth', 1);
        clearvars vp vs rho Vw Is;
        hold on
    end
end
hold on
%% Vli
figure('Color','w'); 
grid on;
scatter3(ISob, VRobs, Den, 20, Vli, ...
    'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat')
xlabel('Shear Impedance'); ylabel('Vp/Vs Ratio'); zlabel('Density');
xlim([3 7]); ylim([1.2 2]); zlim([2.3 2.7]);
view(-12, -43);colormap(jet);colorbar;hold on 
for k=1:length(vli)
    for j=1:length(por)
        for i=1:length(xp)
vsh=1-vli(k)-0.075;
vqu=1-vsh-vli(k)-0.025;
[vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(j), Den_mean, xp(i));
        end
        Vw(:,j,k)=vp(:,j,k)./vs(:,j,k);
        Is(:,j,k)=vs(:,j,k).*rho(:,j,k)./1000;  
        plot3(Is(:,j,k),Vw(:,j,k),rho(:,j,k), 'k','LineWidth', 1);
        hold on
    end
hold on
    for j=1:length(xp)
        for i=1:length(por)
vsh=1-vli(k)-0.075;
vqu=1-vsh-vli(k)-0.025;
vqu(vqu<0)=0.001;
[vp(j,i,k), vs(j,i,k), rho(j,i,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(i), Den_mean, xp(j));
        end
        Vw(j,:,k)=vp(j,:,k)./vs(j,:,k); 
        Is(j,:,k)=vs(j,:,k).*rho(j,:,k)./1000;  
        plot3(Is(j,:,k),Vw(j,:,k),rho(j,:,k), 'k', 'LineWidth',1);
        hold on
    end
end
hold on
for k=1:length(por)
    for j=1:length(vli)
        for i=1:length(xp)
vsh=1-vli(j)-0.075;
vqu=1-vsh-vli(j)-0.025;
vqu(vqu<0)=0.001;
[vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(j), TOC_mean, SW_mean, por(k), Den_mean, xp(i));
        end
        Vw(:,j,k)=vp(:,j,k)./vs(:,j,k);
        Is(:,j,k)=vs(:,j,k).*rho(:,j,k)./1000;   
        plot3(Is(:,j,k),Vw(:,j,k),rho(:,j,k), 'k', 'LineWidth', 1);
        clearvars vp vs rho Vw Is;
        hold on
    end
hold on
    for j=1:length(xp)
        for i=1:length(vli)
vsh=1-vli(i)-0.075;
vqu=1-vsh-vli(i)-0.025;
vqu(vqu<0)=0.001;
[vp(j,i,k), vs(j,i,k), rho(j,i,k)] = modeling(vsh, vqu, vli(i), TOC_mean, SW_mean, por(k), Den_mean, xp(j));
        end
        Vw(j,:,k)=vp(j,:,k)./vs(j,:,k); 
        Is(j,:,k)=vs(j,:,k).*rho(j,:,k)./1000;    
        plot3(Is(j,:,k),Vw(j,:,k),rho(j,:,k), 'k', 'LineWidth', 1);
        clearvars vp vs rho Vw Is;
        hold on
    end
end
hold on
%% XP
figure('Color','w'); 
scatter3(ISob, VRobs, Den, 20, XP, ...
    'Marker', 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat')
xlabel('Shear Impedance'); ylabel('Vp/Vs Ratio'); zlabel('Density');
xlim([3 7]); ylim([1.2 2]); zlim([2.3 2.7]);
grid on;view(55, -7);colormap(jet);colorbar;hold on 
for k=1:length(vli)
    for j=1:length(por)
        for i=1:length(xp)
vsh=1-vli(k)-0.075;
vqu=1-vsh-vli(k)-0.025;
[vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(j), Den_mean, xp(i));
        end
        Vw(:,j,k)=vp(:,j,k)./vs(:,j,k); 
        Is(:,j,k)=vs(:,j,k).*rho(:,j,k)./1000;   
        plot3(Is(:,j,k),Vw(:,j,k),rho(:,j,k), 'k','LineWidth', 1);
        hold on
    end
hold on
    for j=1:length(xp)
        for i=1:length(por)
vsh=1-vli(k)-0.075;
vqu=1-vsh-vli(k)-0.025;
vqu(vqu<0)=0.001;
[vp(j,i,k), vs(j,i,k), rho(j,i,k)] = modeling(vsh, vqu, vli(k), TOC_mean, SW_mean, por(i), Den_mean, xp(j));
        end
        Vw(j,:,k)=vp(j,:,k)./vs(j,:,k);
        Is(j,:,k)=vs(j,:,k).*rho(j,:,k)./1000;   
        plot3(Is(j,:,k),Vw(j,:,k),rho(j,:,k), 'k', 'LineWidth',1);
        hold on
    end
end
hold on
for k=1:length(por)
    for j=1:length(vli)
        for i=1:length(xp)
vsh=1-vli(j)-0.075;
vqu=1-vsh-vli(j)-0.025;
vqu(vqu<0)=0.001;
[vp(i,j,k), vs(i,j,k), rho(i,j,k)] = modeling(vsh, vqu, vli(j), TOC_mean, SW_mean, por(k), Den_mean, xp(i));
        end
        Vw(:,j,k)=vp(:,j,k)./vs(:,j,k); 
        Is(:,j,k)=vs(:,j,k).*rho(:,j,k)./1000;   
        plot3(Is(:,j,k),Vw(:,j,k),rho(:,j,k), 'k', 'LineWidth', 1);
        clearvars vp vs rho Vw Is;
        hold on
    end
hold on
    for j=1:length(xp)
        for i=1:length(vli)
vsh=1-vli(i)-0.075;
vqu=1-vsh-vli(i)-0.025;
vqu(vqu<0)=0.001;
[vp(j,i,k), vs(j,i,k), rho(j,i,k)] = modeling(vsh, vqu, vli(i), TOC_mean, SW_mean, por(k), Den_mean, xp(j));
        end
        Vw(j,:,k)=vp(j,:,k)./vs(j,:,k); 
        Is(j,:,k)=vs(j,:,k).*rho(j,:,k)./1000;   
        plot3(Is(j,:,k),Vw(j,:,k),rho(j,:,k), 'k', 'LineWidth', 1);
        clearvars vp vs rho Vw Is;
        hold on
    end
end
