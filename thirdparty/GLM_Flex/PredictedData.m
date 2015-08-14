function p = PredictedData(y,x)

x = [ones(size(x,1)) x];
b = pinv(x)*y;

p = x*b;