function BlockWrapper(Subj, ord);
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
    ord = [1 2 3];
end

mkdir([Subj '/functional/First_Level_Models/BlockAna_Standard']);
try
    copyfile([Subj '/functional/Standard_Preproc/RunIDs.txt'], [Subj '/functional/First_Level_Models/BlockAna_Standard/']);
catch
    system(['cp ' [Subj '/functional/Standard_Preproc/RunIDs.txt'] ' ' [Subj '/functional/First_Level_Models/BlockAna_Standard/']]);
end
tmp = {'Run' 'Condition' 'Onset' 'Duration'
    1   'Novel'     5       40
    1   'Repeated'  70      40
    1   'Novel'     135     40
    1   'Repeated'  200     40
    2   'Novel'     5       40
    2   'Repeated'  70      40
    2   'Novel'     135     40
    2   'Repeated'  200     40
    3   'Novel'     5       40
    3   'Repeated'  70      40
    3   'Novel'     135     40
    3   'Repeated'  200     40
    4   'Novel'     5       40
    4   'Repeated'  70      40
    4   'Novel'     135     40
    4   'Repeated'  200     40
    5   'Novel'     5       40
    5   'Repeated'  70      40
    5   'Novel'     135     40
    5   'Repeated'  200     40
    6   'Novel'     5       40
    6   'Repeated'  70      40
    6   'Novel'     135     40
    6   'Repeated'  200     40};

WriteDataToText(tmp, [Subj '/functional/First_Level_Models/BlockAna_Standard/RunInfo.csv'], 'w', ',');
SPM8_FirstLevel(Subj,'MixedDes_Block_Params', ord);

