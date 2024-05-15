function stars = significance_text(p,firststarp)
    if ~exist('firststarp','var') || isempty(firststarp)
        firststarp = 0.05;
    end
    if p<=1E-3 && firststarp>=1E-3
        stars='***'; 
    elseif p<=1E-2 && firststarp>=1E-2
        stars='**';
    elseif p<=0.05 && firststarp>=0.05
        stars='*';
    elseif p<=0.2
        stars=['p=' num2str(p,4)];
    elseif isnan(p)
        stars='NaN';
        warning('NaN p value');
    else
        stars='';
    end
end


