function [y] = isrowvector(x)
% [y] = isrowvector(x)
%
% written by ESBM

y = isvector(x) && size(x,1) == 1;