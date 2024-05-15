function [p] = combine_pvalues(pvals,dim)
%combine_pvalues(pvals)
%combine_pvalues(pvals,dim)
%computes a p-value from the product of the input p-values
%  along the dimension 'dim'
%  (if pvals is a vector, defaults to a produce along its length;
%   if pvals is a matrix, defaults to dimension = 1)
%
%e.g. say you have three datasets where condition A has a
% bigger mean than condition B, but only at p=.1 for each case.
% No single dataset is significant. But you don't think that
% A and B are the same, because if they were, it would be very unlikely to 
% get such low p-values three times in a row. You need a way to combine the
% three p-values together into an 'overall' p-value.
% 
% Fisher defined a test statistic for this, 
% -2*ln(p1*p2*p3*...*pn)
% which if all null hypotheses were true, is said to have a chi-square 
% distribution with 2*n degrees of freedom.
%
% See here:
% Fisher, R. A., Statistical Methods for Research Workers, London: Oliver and Boyd, 11th ed., 1950, 
% pages 103-105.
% http://www.haghish.com/resources/materials/Statistical_Methods_for_Research_Workers.pdf
%
% written by ESBM

if nargin < 2
    if isrowvector(pvals)
        dim = 2;
    elseif iscolvector(pvals)
        dim = 1;
    else
        dim = 1;
    end;
end;

p = 1 - chi2cdf( -2*sum(log(pvals),dim), 2*size(pvals,dim));
