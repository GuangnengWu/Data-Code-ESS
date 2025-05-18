function [mupost, sigmapost, Ppost]  = RockPhysicsInversion(mum, sm, G, mdomain, dcond, sigmaerr)

if isrow(mum)
    mum = mum';
end

nv = length(mum);
ns = size(dcond,1);
% Compute prior predicted data mean and covariances
mud = G * mum;
sd  = G * sm * G';
smd = sm * G';
sdm = G * sm;
% Compute posterior covariance of model parameters
sigmapost = sm - smd / (sd + sigmaerr) * sdm;
sigmapost = (sigmapost + sigmapost') / 2;
sigmapost = sigmapost + 1e-6 * eye(size(sigmapost));

mupost = zeros(ns, nv);
Ppost  = zeros(ns, size(mdomain, 1));

for i = 1:ns
    mu_i = mum + smd / (sd + sigmaerr) * (dcond(i,:)' - mud); 
    mupost(i,:) = mu_i'; 
    Ppost(i,:) = mvnpdf(mdomain, mu_i', sigmapost); 
end