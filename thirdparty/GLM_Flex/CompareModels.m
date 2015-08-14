function OUT = CompareModels(mod1,mod2)
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

y1 = sort(mod1.ModelSummary.Y);
y2 = sort(mod2.ModelSummary.Y);

if (numel(y1)~=numel(y2)) || ~all(y1==y2)
    OUT = 'Go to hell!';
    warning('Data between models is not identical!!! Aborting now...');
    return
end

RSS = [mod1.ModelSummary.RSS mod2.ModelSummary.RSS];
df = [mod1.ModelSummary.df(2) mod2.ModelSummary.df(2)];

[trash ord] = sort(df);

RSS = RSS(ord);
df = df(ord);

ndf = diff(df);
SSC = diff(RSS);

F = (diff(RSS)/ndf) / (RSS(1)/df(1));
p = 1-spm_Fcdf(F,ndf, df(1));


OUT.RSS = RSS;
OUT.df = df;
OUT.SSC = SSC;
OUT.F = F;
OUT.ndf = [ndf df(1)];

name = 'Model Difference';
spacer = numel(name)+4;
fprintf(['\n%-' num2str(spacer) 's%s\t%s\t%s\t%s\t%s\n'],'','df','Sum Sq', 'Mean Sq','F-value', 'P-value');
fprintf(['%-' num2str(spacer) 's%0.0f\t%6.3f\t%6.3f\t%6.3f\t%1.6f\n'],name,ndf,SSC,SSC/ndf,F,p);
fprintf(['%-' num2str(spacer) 's%0.0f\t%6.3f\t%6.3f\n'],'Error',df(1),RSS(1),RSS(1)/df(1));


