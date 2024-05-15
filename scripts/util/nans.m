function [n] = nans(varargin)
% [n] = nans(...)
% returns a matrix of NaNs
% nans(...) is the same as nan*ones(...)
%
% written by ESBM

n = nan*ones(varargin{:});
