function bspm_check_motion_print(rpfiles)
% BSPM_CHECK_MOTION Checks motion descriptives form an rp*txt file
% 
% USAGE: bspm_check_motion(rpfiles)
%
% ARGUMENTS
%   rpfiles = realignment parameter files (cell arary)
%

% ------------------------ Copyright (C) 2014 ------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, mfile_showhelp; return; end
if ischar(rpfiles), rpfiles = cellstr(rpfiles); end
for i = 1:length(rpfiles)
    [runpath rpname e] = fileparts(rpfiles{i});
    [rawpath runname{i} e] = fileparts(runpath);
    [subpath rawname e] = fileparts(rawpath);
    [studypath subname{i} e] = fileparts(subpath);
end
sublist = unique(subname);
for i = 1:length(sublist)
    
    idx = cellstrfind(subname,sublist{i});
    crun = runname(idx);
    crp = rpfiles(idx);
    figure('Color','white');
    for r = 1:length(crp)
        
        subplot(1,length(crun),r);
        rp = load(crp{r});
        rp(:,4:6) = rp(:,4:6)*57.3;
        plot(rp, 'LineWidth', .75);
        title(sprintf('%s',regexprep(crun{r},'_',' ')));
        axis([1 length(rp) min(rp(:)) max(rp(:))]);
        set(get(gca,'YLabel'), 'String', 'Movement (mm/degrees)', 'FontSize', 10);
        set(get(gca,'XLabel'), 'String', 'Scan', 'FontSize', 10);
%         legendh = legend({'x' 'y' 'z' 'pitch' 'roll' 'yaw'});
%         set(legendh, 'Orientation', 'vertical', 'Location', 'bestoutside', 'FontSize', 12);
        
    end
    outname = sprintf('%s/MOTION_%s.pdf',studypath,sublist{i});
    
    % resize figure
    dim = [length(crp) .9]*4;
    set(gcf,'units','inches','position',[0 0 dim(1) dim(2)-1]);
    
    % save it
    export_fig(outname);
    close all
    
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
