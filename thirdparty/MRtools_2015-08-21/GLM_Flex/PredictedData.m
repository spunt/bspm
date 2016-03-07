function [p R2] = PredictedData(y,x)

x = [ones(size(x,1)) x];
b = pinv(x)*y;

p = x*b;
R2 = corr(p,y).^2;