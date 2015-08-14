function dv = makedummy(input,difference)
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

if nargin == 1;
    difference = 0;
end
   
if numel(size(input))>2
    error('This function can only be passed a vector');
end

if ~any(size(input)==1)
    error('This function can only be passed a vector');
end

input = input(:);

if isobject(input) && iscolumn(input)
    input = cellstr(input);
end

if iscell(input)
    list= unique(input,'stable');
    dv = zeros(size(input,1),numel(list));
    for ii = 1:numel(list);
        i1 = strmatch(list{ii}, input, 'exact');
        dv(i1,ii) = 1;
    end
end


if isnumeric(input);
    list= unique(input,'stable');
    dv = zeros(size(input,1),numel(list));
    for ii = 1:numel(list);
        i1 = find(input==list(ii));
        dv(i1,ii) = 1;
    end
end

if difference
    tmp = -1*diff(dv,1,2);
    if ~isempty(tmp)
        dv = tmp;
    end
end