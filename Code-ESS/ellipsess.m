function [P, Q] = ellipsess(Km, um, Ki, ui, arf)
% Calculate P and Q values for oblate spheroidal inclusions
% Inputs:
%   Km  - Matrix bulk modulus (from VRH averaging)
%   um  - Matrix shear modulus (from VRH averaging)
%   Ki  - Inclusion bulk modulus (from experiment or empirical data)
%   ui  - Inclusion shear modulus (from experiment or empirical data)
%   arf - Aspect ratio of the ellipsoidal inclusion
% Outputs:
%   P, Q - Parameters used in effective medium theory (Eshelby-based)

if (arf >= 1)
    st = arf ./ (arf.^2 - 1).^1.5 .* (arf .* sqrt(arf.^2 - 1) - acosh(arf));
else
    st = arf ./ (1 - arf.^2).^1.5 .* (acos(arf) - arf .* sqrt(1 - arf.^2));
end

f = arf.^2 ./ (1 - arf.^2) .* (3 .* st - 2);
A = ui ./ um - 1;
B = 1/3 * (Ki ./ Km - ui ./ um);
vm = (3 .* Km - 2 .* um) ./ (2 .* (3 .* Km + um));
R = (1 - 2 .* vm) ./ (2 .* (1 - vm));

F1 = 1 + A .* (1.5 .* (f + st) - R .* (1.5 .* f + 2.5 .* st - 4/3));
F2 = 1 + A .* (1 + 1.5 .* (f + st) - R / 2 .* (3 .* f + 5 .* st)) + ...
     B .* (3 - 4 .* R) + A / 2 .* (A + 3 .* B) .* (3 - 4 .* R) .* ...
     (f + st - R .* (f - st + 2 .* st.^2));
F3 = 1 + A .* (1 - (f + 1.5 .* st) + R .* (f + st));
F4 = 1 + A / 4 .* (f + 3 .* st - R .* (f - st));
F5 = A .* (-f + R .* (f + st - 4/3)) + B .* st .* (3 - 4 .* R);
F6 = 1 + A .* (1 + f - R .* (f + st)) + B .* (1 - st) .* (3 - 4 .* R);
F7 = 2 + A / 4 .* (3 .* f + 9 .* st - R .* (3 .* f + 5 .* st)) + ...
     B .* st .* (3 - 4 .* R);
F8 = A .* (1 - 2 .* R + f / 2 .* (R - 1) + st / 2 .* (5 .* R - 3)) + ...
     B .* (1 - st) .* (3 - 4 .* R);
F9 = A .* ((R - 1) .* f - R .* st) + B .* st .* (3 - 4 .* R);

T1 = 3 .* F1 ./ F2;
T2 = 1/3 .* T1 + 2 ./ F3 + 1 ./ F4 + (F4 .* F5 + F6 .* F7 - F8 .* F9) ./ (F2 .* F4);

P = 1/3 .* T1;
Q = 1/5 .* (T2 - 1/3 .* T1);