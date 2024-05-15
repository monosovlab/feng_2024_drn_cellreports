function val = ternary(bool,true_val,false_val)
% function val = ternary(bool,true_val,false_val)
%   ternary operator

if (bool)
    if iscell(true_val) && length(true_val) == 2 && ((isstring(true_val{2}) && length(true_val{2} == 1)) || ischar(true_val{2}))
        val = field_or_NaN(true_val{1},true_val{2});
    else
        val = true_val;
    end
else
    if iscell(false_val) && length(false_val) == 2 && ((isstring(false_val{2}) && length(false_val{2} == 1)) || ischar(false_val{2}))
        val = field_or_NaN(false_val{1},false_val{2});
    else
        val = false_val;
    end
end