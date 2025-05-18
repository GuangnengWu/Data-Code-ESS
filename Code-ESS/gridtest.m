function [vp_pre_best, vs_pre_best, xp_best] = gridtest(VP, VS, Vsh, Vqu, Vli, TOC, Den, SW, POR)
% Grid search method for inverting the pore aspect ratio (xp) 
% by minimizing the error between modeled and observed VP

tol = inf;  % Initial error set to infinity
xp_best = NaN;
vp_pre_best = NaN;
vs_pre_best = NaN;

for xp = 0.03:0.01:0.99
    [vp_pre, vs_pre, den] = modeling(Vsh, Vqu, Vli, TOC, SW, POR, Den, xp);

    % Mean squared error between predicted and observed VP
    err = mean((vp_pre(:) - VP(:)).^2);
    
    % Alternative combined error using both VP and VS (optional)
    % err = abs(vp_pre - VP) ./ VP + abs(vs_pre - VS) ./ VS;

    if err < tol
        tol = err;
        vp_pre_best = vp_pre;
        vs_pre_best = vs_pre;
        xp_best = xp;
    end
end

end