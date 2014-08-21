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

cols = {'Max Trans' 'Max Rot' 'Max Diff'};

for i = 1:length(rpfiles)
        
    [path name ext] = fileparts(rpfiles{i});
    [path name ext] = fileparts(path);
    run{i} = name;
    [path name ext] = fileparts(path);
    [path name ext] = fileparts(path);
    sub{i} = name;
    rp = load(rpfiles{i});
    rp(:,4:6) = rp(:,4:6)*57.3;
    rpdiff = diff(rp);
    movement = max(rp) - min(rp);
    maxtrans(i) = max(movement(1:3));
    maxrot(i) = max(movement(4:6));
    maxdiff(i) = max(max(rpdiff));
    
end
data = [maxtrans' maxrot' maxdiff'];
out.data = data; 
out.cols = cols; 
out.subs = sub; 
gdata = [nanmean(data); nanstd(data); max(data)];
data = [data; gdata];
gnames = {'MEAN' 'SD' 'MAX'};
sub = [sub gnames];
disptable(data,cols,sub,'%2.2f');

if writereport
    [d t] = bob_timestamp;
    outname = sprintf('motion_report_%s_%s.xls',d,t);
    allcell = [[{''} cols]; [sub' num2cell(data)]];
    xlwrite(outname, allcell);
end



    
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
