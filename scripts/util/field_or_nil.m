function y = field_or_nil(struct, field, nilval)
% y = field_or_NaN(struct, field, nilval)
%   if struct has field, return val, otherwise just return nilval (default
%   is NaN) 

    if isfield(struct, field)
        y = struct.(field);
    else
        if ~exist('nilval','var')
            nilval = NaN;
        end
        y = nilval;
    end
end
