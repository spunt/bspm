function [pc,perc,available_number] = bramila_princomp(ts,requested_number)

% making sure that ts has no NaNs or Inf since svd does note like those
ts(find(isnan(ts)))=0;
ts(find(isinf(ts)))=0;

if nargin>1 % internal PCA (signs might be flipped)
    warning('off', 'stats:pca:ColRankDefX');
    [~,score,~,~,explained]= pca(ts,'NumComponents',requested_number,'Centered',true);
    available_number=size(score,2);
    if available_number<requested_number
        warning('Number of available components is less than requested (%i<%i)',available_number,requested_number);
    end
    perc=cumsum(explained);
    perc = perc(1:available_number);
    pc=score(:,1:available_number);
else % direct SVD (for first PC only, correct sign)
    [m n]   = size(ts);
    if m > n
        [v,s,v] = svd(ts'*ts); %#ok<*ASGLU>
        s = diag(s);
        v = v(:,1);
        u = ts*v/sqrt(s(1));
    else
        [u,s,u] = svd(ts*ts');
        s       = diag(s);
        u       = u(:,1);
        v       = ts'*u/sqrt(s(1));
    end
    d  = sign(sum(v));
    u  = u*d;
    v  = v*d;
    pc = u*sqrt(s(1)/n);
    perc=100*cumsum(s)/sum(s);
end

end
