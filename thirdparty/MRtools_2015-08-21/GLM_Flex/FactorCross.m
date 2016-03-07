function X = FactorCross(varargin)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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

go = 1; 

input = [];
N = [];
if numel(varargin) == 1
    input = varargin{1};
    for ii = 1:numel(input)
        N(ii) = size(input{ii},2);
    end
else
    for ii = 1:numel(varargin)
        input{ii} = varargin{ii};
        N(ii) = size(input{ii},2);
    end
end

X = zeros(size(input{1},1),prod(N));

counts = ones(size(N));
allCounts = []; c = 0;
while go
    
    st = ones(size(input{1},1),1);
    for ii = 1:numel(N);
        st = st.*input{ii}(:,counts(ii));
    end
   
    
    %disp(counts)
    for ii = numel(counts):-1:0;
        if ii == 0;
            go=0;
            break
        end
        if counts(ii)==N(ii);
            counts(ii) = 1;
        else
            counts(ii) = counts(ii)+1;
            break;
        end
    end
   
    
    c = c+1;
    allCounts(c,:) = counts;
    X(:,c) = st;
    
end
    

