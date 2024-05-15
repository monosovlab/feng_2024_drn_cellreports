function escaped_strs = escape(strs, char, escapechar)
% escaped_strs = escape(strs, char)
% 	- add a backslash before each character <char> in each string/char array in <strs>
    if ~exist('char','var')
        char = '_';
    end
    if ~exist('escapechar','var')
        escapechar = '\\';
    end
	replace = string([escapechar,char]);
	if iscell(strs)
		escaped_strs = cellfun(@(str) regexprep(str,char,replace), strs,'UniformOutput',false);
	else
		escaped_strs = regexprep(strs,char,replace);
	end
end
