function [Y ch1 ch2] = QuickROI(mni,ss, fn);
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
list=fn;

clear Y;
for ii = 1:length(list)
    D = getROI(list{ii}, ss, mni);
%     Y2(ii) =   D.svd;
    Y(ii) =  D.mean;
    ch1(ii) = D.nNaN(1);
    ch2(ii) = D.nNaN(2);
end
% keyboard;