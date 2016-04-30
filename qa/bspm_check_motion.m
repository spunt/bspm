function out = bspm_check_motion(rpfiles, writereport)
% BSPM_CHECK_MOTION Checks motion descriptives form an rp*txt file
% 
% USAGE: out = bspm_check_motion(rpfiles, writereport)
%
% ARGUMENTS
%   rpfiles = realignment parameter files (cell arary)
%   writereport = option to write out report
%

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, writereport = 0; end
if nargin<1, display('USAGE: out = bspm_check_motion(rpfiles, writereport)'); return; end
if ischar(rpfiles), rpfiles = cellstr(rpfiles); end
nrp     = length(rpfiles);
data    = NaN(nrp, 4);
subcell = cell(nrp, 3);
for i = 1:nrp
    
    [path name ext] = fileparts(rpfiles{i});
    subcell{i,3} = name; 
    [path name ext] = fileparts(path);
    subcell{i,2} = name;
    [path name ext] = fileparts(path);
    [path name ext] = fileparts(path);
    subcell{i,1} = name;
    cdata = load(rpfiles{i});
    rp = cdata(:,1:6);
%     rp(:,4:6) = rp(:,4:6)*57.3;
    rpdiff = diff(rp);
    movement = max(rp) - min(rp);
    data(i,1) = max(movement(1:3));
    data(i,2) = max(movement(4:6));
    data(i,3) = max(max(rpdiff));
    if size(cdata, 2) > 6
        dumdim = size(cdata(:,7:end));
        data(i,4) = 100*(dumdim(2)/dumdim(1));
    end
end
cols = {'Subject' 'Run' 'Filename' 'Max Trans' 'Max Rot' 'Max Diff'};
if all(isnan(data(:,4)))
    data(:,4) = []; 
else
    cols = [cols {'% Bad'}]; 
end
allcell = [cols; [subcell num2cell(data)]];
if writereport
    outname = sprintf('motion_report_%s.xlsx', bspm_timestamp);
    xlwrite(outname, allcell);
end
out.data = data; 
out.cols = cols; 
out.subs = subcell;
out.allcell = allcell; 
cols = [cols(4:end) strcat('Z-', cols(4:end))]
disptable([data oneoutzscore(data)],cols,subcell(:,1),'%2.2f');



    
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
