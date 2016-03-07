function [rr ind h] = MakeConnMatrix(fn,dim,greyThresh);
%%% Written by Aaron Schultz (aschultz@martinos.org)
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
%%% GNU General Public License for more details.

if nargin<2
    dim = [4 4 4];
end

if nargin<3
    greyThresh = .3;
end


disp('Reading and Resizing Data');

[pth filename extension] = fileparts(fn);
h = spm_vol([fn ',1']);
[data mat] = SliceAndDice(fn,h,dim,h,[1 NaN],[]);
hm = spm_vol(which('grey.nii'));
[mask mat] = SliceAndDice(hm,h,dim,h,[1 NaN],[]);
ind = find(mask>greyThresh);

h.dim = size(mask); h.mat = mat; h.dt = [16 0];

M = reshapeWholeBrain(size(data),data);

M = zscore(M(:,ind));
rr = zeros(numel(ind),numel(ind),'single');


% q = ReadInFile('/proc/cpuinfo','\t');
% ncores = numel(strmatch('processor',q(:,1)))-1;
% if ncores>8
%     ncores = 6;
% end
%     
% disp('Creating the correlation matrix');
% matlabpool(ncores)
% parfor ii = 1:length(ind)
% % for ii = 1:length(ind)
% %     if mod(ii,1000)==0
% %        disp([ii numel(ind)]); 
% %     end
% %     tmp = '000000000000';
% %     lab = num2str(ind(ii));
% %     tmp(end-length(lab)+1:end)=lab;
% %     eval(['b' tmp ' = single(pinv([M(:,ind(ii))])*M);'])
% %     save('HubDat.mat',['b' tmp],'-append')
%     rr(ii,:) = single(pinv(M(:,ind(ii)))*M(:,ind));
% end
% matlabpool close

% disp('Computing Distances');
% x = [];[x(:,1) x(:,2) x(:,3)] = ind2sub(h.dim,ind);
% D = DistMatF(x,0);

disp('Creating the correlation matrix');
cols = numel(ind);
rows = size(M,1);
blocks=300;
CC = (1:blocks:cols)'; CC(:,2) = [CC(1:end-1)+blocks-1; cols];
% if matlabpool('size') > 0
%     parfor ii=1:size(CC,1);
%         %ti = (((ii-1)*300)+1):(ii*300);
%         %rr(ti,:) = single((M(:,ti)' * M) / (rows-1));
%         rr(CC(ii,1):CC(ii,2),:) = single((M(:,CC(ii,1):CC(ii,2))' * M) / (rows-1));
%     end
% else
    for ii=1:size(CC,1);
        rr(CC(ii,1):CC(ii,2),:) = single((M(:,CC(ii,1):CC(ii,2))' * M) / (rows-1));
    end
% end

% if nargin == 1
%     disp('Saving Data');
%     save([filename '_ConMat.mat'],'rr','ind','h', '-v7.3');
% end