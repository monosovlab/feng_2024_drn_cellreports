function [y] = iscolvector(x)
% [y] = iscolvector(x)
%
% written by ESBM

y = isvector(x) && size(x,2) == 1;