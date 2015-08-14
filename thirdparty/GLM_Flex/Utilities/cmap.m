function [cols cm cc] = cmap(X, lims, cm)
%%% Written by Aaron P. Schultz - aschultz@martinos.org
%%%
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.X = X(:);
lims = sort(lims);
nBins = 256;

if ischar(cm)
    %eval(['cm = ' cm '(' num2str(nBins) ');']);
    eval(['cm = colmap(''' cm ''',' num2str(nBins) ');']);
end

X(find(X<lims(1)))=lims(1);
X(find(X>lims(2)))=lims(2);

% bins = (lims(1):diff(lims)/(nBins):lims(2)); 
% % bins(1) = bins(1)+.01; bins(end) = bins(end)-.01
% % bins = [-inf bins inf];
% bb = zeros(size(X));
% for ii = 1:numel(bins)-1;
%     i1 = find(X>=bins(ii) & X<=bins(ii+1));
%     bb(i1) = ii;
% end
% [min(bb) max(bb)]


cc = [X(:); lims(:)];
cc = cc/(diff(lims));
cc = (cc*(nBins));
dd = cc;
cc = cc+(nBins-max(cc));
cc = floor(cc(1:end-2))+1;
cc(cc>nBins)=nBins;
cc(cc<1)=1;
% [min(cc) max(cc)]
% figure(30); plot(bb,cc,'.'); shg

% try
% cols = cm(cc,:);
cols = nan(numel(cc),3);
cols(~isnan(cc),:) = cm(cc(~isnan(cc)),:);
% catch
%     keyboard;
% end
end
function out = cmap_upsample(in, N)
    num = size(in,1);
    ind = repmat(1:num, ceil(N/num), 1);
    rem = numel(ind) - N; 
    if rem, ind(end,end-rem+1:end) = NaN; end
    ind = ind(:); ind(isnan(ind)) = [];
    out = in(ind(:),:);
end

