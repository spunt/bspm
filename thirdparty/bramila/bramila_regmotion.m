

function [voldata r2] = bramila_regmotion(cfg)
ts=cfg.ts;
voldata=cfg.vol;
% function [b r2] = fmriregress2(ts, voldata)
%
% TS = time series
% VOLDATA = volumedata
% B = 4-D matrix containing the regression coefficients
% R2 = 3-D matrix containing the R^2

ts=zscore(ts);
% n_ts=[ones(size(ts,1),1) ts];
s=size(voldata);
% b=zeros(size(n_ts,2),s(1),s(2),s(3));
% res=zeros(size(voldata));
r2=zeros(s(1),s(2),s(3));
%h = waitbar(0,'Hang in there .........');
for k=1:s(1)
fprintf([num2str(k) '..']);
%    waitbar(k/s(1),h)
    for m=1:s(2)
        for n=1:s(3)
            y=squeeze(voldata(k,m,n,:));
	    my=mean(y);
	    y=y-my;
            [bb int R]=regress(y,ts);
%             b(:,k,m,n)=bb;
            voldata(k,m,n,:)=R+my;
            cc=corrcoef([y ts*bb]);
            r2(k,m,n)=cc(1,2)^2;
        end
    end
end
r2(find(isnan(r2)))=0;
%close(h)


