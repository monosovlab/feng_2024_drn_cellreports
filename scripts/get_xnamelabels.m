function xnamelabels = get_xnamelabels(xname)
    xnamelabels = xname;
    xnamelabels = regexprep(xnamelabels,'Rgt','rgt');
    xnamelabels = regexprep(xnamelabels,'R','E[R]');
    xnamelabels = regexprep(xnamelabels,'O','T[R]');
    xnamelabels = regexprep(xnamelabels,'S','Unc[R]');
    xnamelabels = regexprep(xnamelabels,'I','Info');
    xnamelabels = regexprep(xnamelabels,'rgt','Rgt');
end