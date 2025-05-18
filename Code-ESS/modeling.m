function [vp_pre, vs_pre, den] = modeling(Vsh, Vqu, Vli, TOC, SW, POR, Den, xp)

%--------------------------------------------------------------------------
% Moduli of matrix minerals and fluids
% pp = [1 k_quartz  2 k_calcite  3 k_shale  4 k_kerogen  
%       5 u_quartz  6 u_calcite  7 u_shale 8 u_kerogen  
%       9 k_water  10 k_oil]
pp = [37     77     25    2.9    44     32     9     2.7    2.25   0.57 ];

% Densities of matrix minerals and fluids
% den_init = [ quartz  calcite  shale  kerogen  water  oil ]
den_init = [ 2.65     2.71     2.55    1.4      1.04   0.7 ];
%--------------------------------------------------------------------------

% First mixing stage: brittle minerals (quartz and calcite)
V1 = Vqu + Vli;
Vqu1 = Vqu ./ V1;
Vli1 = 1 - Vqu1;
[k0, u0] = bound(0, [Vqu1, Vli1], [pp(1), pp(2)], [pp(5), pp(6)]); % VRH average

% Second mixing stage: ductile minerals (shale and TOC)
V2 = Vsh + TOC;
Vsh1 = Vsh ./ V2;
TOC1 = 1 - Vsh1;
[k1, u1] = kuster(pp(3), pp(7), pp(4), pp(8), 0.1, TOC1); % Kuster-Toksöz

% Normalize volume fractions of all minerals
Vsh2 = Vsh ./ (Vsh + Vli + Vqu + TOC);
Vli2 = Vli ./ (Vsh + Vli + Vqu + TOC);
Vqu2 = Vqu ./ (Vsh + Vli + Vqu + TOC);
TOC2 = TOC ./ (Vsh + Vli + Vqu + TOC);

% Densities of brittle and ductile mixtures
den0 = Vli1 .* den_init(2) + Vqu1 .* den_init(1); % Brittle mineral density
den1 = Vsh1 .* den_init(3) + TOC1 .* den_init(4); % Ductile mixture density

% P-wave and S-wave velocities of brittle and ductile parts (in m/s)
vp_0 = sqrt((k0 + 4/3 .* u0) .* 1e6 ./ den0);
vs_0 = sqrt(u0 .* 1e6 ./ den0);
vp_1 = sqrt((k1 + 4/3 .* u1) .* 1e6 ./ den1);
vs_1 = sqrt(u1 .* 1e6 ./ den1);

% Backus average to get anisotropic stiffness constants
[c11_bk, c12_bk, c13_bk, c33_bk, c44_bk, c, den_0] = bkusc( ...
    [Vli2 + Vqu2, Vsh2 + TOC2], [vp_0, vp_1], [vs_0, vs_1], [den0, den1]);

%--------------------------------------------------------------------------
% Porosity normalization and effective saturated rock density
% denF=1.05;
% pt=(Den-den_0)./(denF-den_0);
% POR=pt;
POR1 = POR ./ (Vsh + Vli + Vqu + TOC + POR);
den = (1 - POR1) .* den_0 + POR1 .* SW .* den_init(5) + POR1 .* (1 - SW) .* den_init(6);

%--------------------------------------------------------------------------
% Stiffness of dry rock using SCA-DEM model
c1 = zeros(6, 6); % Initial dry rock stiffness
[c11_sca, c33_sca, c44_sca, c12_sca, c13_sca] = sca(c, c1, xp, 1 - POR1, POR1);

%--------------------------------------------------------------------------
% Fluid bulk modulus using Wood’s formula
kf = 1 ./ (SW ./ pp(9) + (1 - SW) ./ pp(10));

%--------------------------------------------------------------------------
% Brown and Korringa theory (anisotropic Gassmann fluid substitution)
% Mineral     - sso (compliance of mineral frame)
% Dry rock    - ssd (compliance of dry rock)
% Saturated   - ssw (compliance of saturated rock)
sso = [1 ./ c11_bk, 1 ./ c12_bk, 1 ./ c13_bk, 1 ./ c33_bk, 1 ./ c44_bk];     % Mineral frame compliance
ssd = [1 ./ c11_sca, 1 ./ c12_sca, 1 ./ c13_sca, 1 ./ c33_sca, 1 ./ c44_sca]; % Dry rock compliance

ssw = bkti(POR1, 1 ./ kf, sso, ssd); % Fluid substitution

% Final predicted velocities (m/s)
vp_pre = sqrt(1 ./ (ssw(4) .* den));
vs_pre = sqrt(1 ./ (ssw(5) .* den));
end