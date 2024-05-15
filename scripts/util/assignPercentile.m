% its actually quantile
function Y = assignPercentile(X)
    Y = (floor(tiedrank(X))-0.5)/length(X);
end

%% test
% tol = 1e-8;
% yfrand = rand(size(1:100));
% yfrandrep = yfrand([1:10 1:10 21:end]);
% all(abs(quantile(yfrandrep,assignPercentile(yfrandrep))-yfrandrep)<tol);