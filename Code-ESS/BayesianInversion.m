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

% Forward modeling to generate rock physics templates
for k = 1:length(vli)
    for j = 1:length(por)
        for i = 1:length(xp)
            vsh = max(0, 1 - vli(k) - 0.075);
            vqu = max(0, 1 - vsh - vli(k) - 0.025);
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

% Flatten synthetic elastic parameters to 1D vectors
x = Is;   xx = x(:);
y = Vw;   yy = y(:);
z = rho;  zz = z(:);

% Flatten modeled physical properties
aa = lime;     a = aa(:);
bb = porosity; b = bb(:);
cc = poar;     c = cc(:);

% Construct interpolants
F1 = scatteredInterpolant(xx, yy, zz ,a, 'linear', 'linear');
F2 = scatteredInterpolant(xx, yy, zz ,b, 'linear', 'linear');
F3 = scatteredInterpolant(xx, yy, zz ,c, 'linear', 'linear');

% Predict physical properties from interpolants
vlipre = F1(ISob , VRobs , Den); vlipre(vlipre < 0.001) = 0.001;
porpre = F2(ISob , VRobs , Den); porpre(porpre < 0.001) = 0.001;
xppre  = F3(ISob , VRobs , Den); xppre(xppre < 0.001) = 0.001;  xppre(xppre > 1) = 1;

%data-driven constraints
mtrain = [POR, Vli, XP];
nv = size(mtrain, 2);
dtrain = [vlipre, porpre, xppre];
nd = size(dtrain, 2);

% Define physical property domains
PORdomain = 0:0.005:0.4;
Vlidomain  = 0:0.005:1;
xpdomain  = 0:0.01:1;
[P, V, S] = ndgrid(PORdomain, Vlidomain, xpdomain);
mdomain = [P(:), V(:), S(:)];

dcond = dtrain; 
ns = size(dcond,1);

R = zeros(nd, nv+1);
X = [mtrain, ones(ns,1)];
R(1,:) = regress(vlipre, X);
R(2,:) = regress(porpre, X);
R(3,:) = regress(xppre, X);

% Assume measurement noise covariance
sigmaerr = 1e-2 * eye(nd);

% Prior mean and covariance
mum = mean(mtrain);
sm  = cov(mtrain);

G = R(:, 1:nv);
datacond = dcond - R(:, end)';

% Run Bayesian rock physics inversion
[mupost, sigmapost, Ppost] = RockPhysicsInversion(mum, sm, G, mdomain, datacond, sigmaerr);

% Posterior mean and confidence intervals
PORpost = mupost(:,1);
Vlipost  = mupost(:,2);
xppost  = mupost(:,3);xppost(xppost<0.001) = 0.001;

PORlp = PORpost - 1.96 * sqrt(sigmapost(1,1));
Vlilp  = Vlipost  - 1.96 * sqrt(sigmapost(2,2));
xplp  = xppost  - 1.96 * sqrt(sigmapost(3,3));

PORup = PORpost + 1.96 * sqrt(sigmapost(1,1));
Vliup  = Vlipost  + 1.96 * sqrt(sigmapost(2,2));
xpup  = xppost  + 1.96 * sqrt(sigmapost(3,3));

% Clamp confidence bounds
PORlp(PORlp<0) = 0;       PORup(PORup>0.4) = 0.4;
Vlilp(Vlilp<0)   = 0;       Vliup(Vliup>1)     = 1;
xplp(xplp<0)   = 0;       xpup(xpup>1)     = 1;

% Compute marginal posterior distributions
PpostPOR = zeros(ns, length(PORdomain));
PpostVli  = zeros(ns, length(Vlidomain));
Ppostxp  = zeros(ns, length(xpdomain));

for i = 1:ns
    Ppostjoint = reshape(Ppost(i,:), length(PORdomain), length(Vlidomain), length(xpdomain));
    PpostPOR(i,:) = sum(sum(Ppostjoint,3),2);
    PpostVli(i,:)  = sum(sum(Ppostjoint,3),1);
    Ppostxp(i,:)  = sum(sum(Ppostjoint,2),1);

    % Normalize
    PpostPOR(i,:) = PpostPOR(i,:) / sum(PpostPOR(i,:));
    PpostVli(i,:)  = PpostVli(i,:)  / sum(PpostVli(i,:));
    Ppostxp(i,:)  = Ppostxp(i,:)  / sum(Ppostxp(i,:));
end

% Plot results
figure('Color','W')
depth = 1:900;

set(gcf, 'Position', [100, 100, 1200, 500])
subplot(1,3,1)
fill([PORlp; flipud(PORup)], [depth'; flipud(depth')], [1 0.8 0.8], 'EdgeColor','none'); hold on;
plot(POR, depth, 'k', 'LineWidth', 2); 
plot(PORpost, depth, 'r--', 'LineWidth', 2); 
xlabel('\phi (v/v)'); ylabel('Depth (m)');xlim([0 0.4]);
legend('95% CI','True','Inverted','Location','best');

subplot(1,3,2)
fill([Vlilp; flipud(Vliup)], [depth'; flipud(depth')], [1 0.8 0.8], 'EdgeColor','none'); hold on;
plot(Vli, depth, 'k', 'LineWidth', 2); 
plot(Vlipost, depth, 'r--', 'LineWidth', 2); 
xlabel('SH (v/v)'); xlim([0 1]);
legend('95% CI','True','Inverted','Location','best');

subplot(1,3,3)
fill([xplp; flipud(xpup)], [depth'; flipud(depth')], [1 0.8 0.8], 'EdgeColor','none'); hold on;
plot(XP, depth, 'k', 'LineWidth', 2); 
plot(xppost, depth, 'r--', 'LineWidth', 2); 
xlabel('xp (v/v)'); xlim([0 1]);
legend('95% CI','True','Inverted','Location','best');