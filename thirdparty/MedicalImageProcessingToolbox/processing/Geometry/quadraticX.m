function out = quadraticX(x)
% x must be a N x K vector, where K is the dimensions and N the number of
% points

N = size(x,1);
K = size(x,2);

out(:,1)= ones(N,1);
out(:,2:K+1)=x;

combinations = nchoosek(1:K,2);
for i=1:size(combinations,1)
    out(:,K+1+i)=x(:,combinations(i,1)).*x(:,combinations(i,2));
end
out(:,(K+1+size(combinations,1)+1):(2*K+1+size(combinations,1)))=x.^2;

end
