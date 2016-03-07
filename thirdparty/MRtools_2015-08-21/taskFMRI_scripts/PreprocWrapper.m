function PreprocWrapper(Subj,ord)
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

if nargin == 1;
    ord = [1 2 3 4];
end

[a b] = dir_wfp([Subj '/functional/Orig/*.nii']);

if ~(exist([Subj '/functional/Standard_Preproc/RunIDs.txt'])>0)
    b = char(b);
    b = b(:,1:5);
    b = cellstr(b);
    disp(b);
    
    mkdir([Subj '/functional/Standard_Preproc/']);
    WriteDataToText(b, [Subj '/functional/Standard_Preproc/RunIDs.txt'], 'w', '\t',0);
end

if numel(b)~=6;
    warning('We do not appear to have 6 good runs. Check the RunIDs file and make sure the runs are specified correctly.');
end

SPM8_Preproc(Subj,'Preproc_Params',ord);
