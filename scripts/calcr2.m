function [efit,r2] = calcr2(efit)
    sse_total = sum((efit.y - mean(efit.y)).^2);
    sse_resid = sum(efit.stats.resid.^2);
    r2 = 1 - (sse_resid ./ sse_total);
    efit.r2 = r2;
end