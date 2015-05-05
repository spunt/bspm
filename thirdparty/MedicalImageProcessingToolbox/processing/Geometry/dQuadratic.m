function dA = dQuadratic(A)
% if the space is of dimension K, then the size of A is 1 + 2*K + (k 2)=3*K
K=2;
while(true)
    if numel(A)==factorial(K)/(factorial(2)*factorial(K-2))+2*K+1
        break;
    end
    K=K+1;
end


dA = diag(2*A(end-K+1:end));
for i=1:K
    for j=i+1:K
        dA(i,j)=A(i+j+1);
        dA(j,i)=A(i+j+1);
    end
end

end
