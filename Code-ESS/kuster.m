function [k,u]=kuster(km,um,ki,ui,a,xi)
% Implementation of the Kuster-Toksoz model
% [k,u] = kuster(km, um, ki, ui, a, xi)
% km, um ----- Matrix (background medium) bulk and shear moduli
% ki, ui ----- Inclusion (pore/fluid/mineral) bulk and shear moduli
% xi --------- Volume fraction of inclusions
% a ---------- Aspect ratio of inclusions
% k , u ------ Output dry rock bulk modulus and shear modulus
[P,Q]=ellipsess(km,um,ki,ui,a);
sm=um./6.*(9.*km+8.*um)./(km+2.*um);
k=(4.*km.*um + 3.*km.^2 + 4.*P.*ki.*um.*xi - 4.*P.*km.*um.*xi)./(3.*km + 4.*um - 3.*P.*ki.*xi + 3.*P.*km.*xi);
u=(sm.*um + um.^2 + Q.*sm.*ui.*xi - Q.*sm.*um.*xi)./(sm + um - Q.*ui.*xi + Q.*um.*xi);