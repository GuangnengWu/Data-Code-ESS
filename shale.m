clear;
clc;
tic;

% poolobj=gcp('nocreate'); % 检查当前是否已有并行计算池
% if isempty(poolobj)
%     poolsize=4; % 设置计算池的大小
%     parpool(poolsize); % 创建并行计算池
% end

%井数据-------------------------------------------------------------------------------------------
load('niu55_1.mat')
load('xp2.mat')
load('VPPRE.mat')
load('TOC1.mat')
load('xp2test.mat')
load('yuce.mat')
Depth=niu55_1(833:1393,1);  
CNL=niu55_1(833:1393,2);      
Den=niu55_1(833:1393,3);     
Den(isnan(Den))=2.52;
GR=niu55_1(833:1393,4);     
Vli=niu55_1(833:1393,5)./100; 
Vli(Vli<0)=0.005;
OTHR=niu55_1(833:1393,6);     
P_wave=niu55_1(833:1393,7);   
VP=niu55_1(833:1393,8);       
POR=niu55_1(833:1393,10)./100;    
S_wave=niu55_1(833:1393,14);  
VS=niu55_1(833:1393,15);      
Vqu=niu55_1(833:1393,16)./100;
Vqu(Vqu<0)=0.005;
Vsh=niu55_1(833:1393,17)./100;
SW=niu55_1(833:1393,19)./100; 
TOC=TOC1(833:1393,1)./100;

%-------------------------------------------------------------------------------------------
for i=1:length(VP)
    i
[vp_pre11(i,1),vs_pre11(i,1),xp2(i,1),den(i,1),vp_wyllie(i,1),pt(i,1),y1(i,1),BI2(i,1)]=grid_alf(VP(i,1),VS(i,1),Vsh(i,1),Vqu(i,1),Vli(i,1),TOC(i,1),Den(i,1),SW(i,1),POR(i,1));

end

%-------------------------------------------------------------------------------------------
vp_pre11=vp_pre11./1000;
vs_pre11=vs_pre11./1000;
VP=VP./1000;
VS=VS./1000;

figure('Color','W')
subplot(171)
plot(smooth(vp_pre11),Depth,'k','LineWidth',1);hold on
plot(smooth(VP),Depth,'b','LineWidth',1);hold on
xlim([2.5 4.5]),ylim([3530 3600]),set(gca,'YDir','reverse');xlabel('VP（km/s）')
subplot(172)
plot(smooth(vs_pre11),Depth,'k','LineWidth',1);hold on
plot(smooth(VS),Depth,'b','LineWidth',1);hold on
xlim([1.5 2.5]),ylim([3530 3600]),set(gca,'YDir','reverse');xlabel('VS（km/s）')
subplot(173)
plot(smooth(xp2),Depth,'k','LineWidth',1);hold on
xlim([0 1]),ylim([3530 3600]),set(gca,'YDir','reverse');xlabel('xp'); 

figure('Color','W')
subplot(111)
scatter(den,pt,'k','filled');hold on
p = polyfit(den, pt,2);
x_fit = linspace(min(den), max(den), 3000);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-r', 'LineWidth', 2);
y_pred = polyval(p, pt);
R = corrcoef(den, y_pred); 
R_squared = R(1, 2)^2;   
fit_equation = sprintf('y = %.2fx^2 + %.2fx + %.2f', p(1), p(2), p(3));
text(min(x_fit), max(y_fit), {fit_equation, sprintf('R^2 = %.2f', R_squared)}, ...
    'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k');

figure('Color','W')
subplot(111)
scatter(xp2,Vsh,'k','filled');hold on
xlim([0 1]),ylim([0 1]);xlabel('xp'); ylabel('Vsh');
figure('Color','W')
subplot(111)
scatter(xp2,POR,'k','filled');hold on
xlim([0 1]),ylim([0 0.15]);xlabel('xp');ylabel('POR');
figure('Color','W')
subplot(111)
scatter(xp2,VP,'k','filled');hold on
xlim([0 1]),ylim([2500 4500]);xlabel('xp');ylabel('VP');

%-------------------------------------------------------------------------------------------

