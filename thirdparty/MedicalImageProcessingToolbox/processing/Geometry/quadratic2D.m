function out = quadratic2D(A,x)
% x must be a N x 2 vector
% A must be a column vector with the values for a0,a1,...a5 so that
% out = a0 +a1*x +a2*y +a3*x*y a4*x^2+a5*y^2

X = [ones(size(x,1),1) x(:,1) x(:,2) x(:,1).*x(:,2) x(:,1).^2 x(:,2).^2 ];

out = X*A;


end
