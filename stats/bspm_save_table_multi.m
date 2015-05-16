function bspm_save_table_multi(allim,intensity,separation,customextentthresh)
% BSPM_SAVE_TABLE_MULTI
%
%   USAGE: bspm_save_table_multi(allim,intensity,separation,mask)
%
%   ARGUMENTS
%       allim:      filename of input statistic image
%       intensity:  intensity threshold to use
%       separation: minimum peak separation (in mm)
%       customextentthresh: by default, will choose fwe cluster corrected
%       extent to use a custom extent threshold, set a value here
%       
%      
% Created April 8, 2013 - Bob Spunt
% Uses code authored by:
% Drs. Donald McLaren & Aaron Schultz (peak_nii.m) 

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<4, customextentthresh = []; end
if nargin<3, error('USAGE: bspm_save_table_multi(allim,intensity,separation,customextentthresh)'); return; end

%% Start it Off %%
headers1 = {'Analysis Name' '' '' '' '' 'MNI Coordinates' '' ''};
headers2 = {'' 'Region Name' 'L/R' 'Extent' 't-value' 'x' 'y' 'z'};
combinedcell = [headers1; headers2];
emptyrow = cell(1,length(headers1));

for i = 1:length(allim)
    
    image = allim{i};
    
    %% get cluster level correction extent %%
    if isempty(customextentthresh)
        cluster = bspm_cluster_correct(image);
    else
        cluster = customextentthresh; 
    end

    %% default structure input to peak_nii %%
    peaknii.thresh = intensity;
    peaknii.cluster = cluster;
    peaknii.out = '';
    peaknii.sign = 'pos';
    peaknii.type = 'T';
    peaknii.voxlimit = [];
    peaknii.separation = separation;
    peaknii.SPM = 1;
    peaknii.conn = [];
    peaknii.mask = [];
    peaknii.df1 = [];
    peaknii.df2 = [];
    peaknii.nearest = 1;
    peaknii.label = 'aal_MNI_V4';
%     peaknii.label = 'HarvardOxford_cortex';
    
    %% run peak_nii %%
    [voxels] = peak_nii(image,peaknii);

    %% clean up data %%
    if ~isempty(voxels)
        labels = voxels{2};
        labels = regexprep(labels,',','');
        peaks = voxels{1};
        peaks(:,6:7) = [];
        side = labels;
        side(:) = {'L'};
        side(peaks(:,3) > 0) = {'R'};
        side(peaks(:,3)==0) = {'L/R'};
        pad = cell(size(side));
        tmp = [pad labels side num2cell(peaks)];
%         datacell = cell(size(tmp,1),length(emptyrow));
        datacell = tmp;
    else
        datacell = emptyrow; 
        datacell{2} = 'No suprathreshold clusters';
    end
    
    %% get names and save %%
    
    [impath imname e] = fileparts(image);
    [p imfolder e] = fileparts(impath);
    acell = emptyrow;
    acell{1} = [imfolder '__' imname];
    
    %% Combined %%
    combinedcell = [combinedcell; emptyrow; acell; datacell];

    %% Cleanup %%
    clear voxels peaks labels
    tmpfile = files([impath filesep '*_peaks_*_clusters.nii']);
    if ~isempty(tmpfile), delete(tmpfile{:}); end
    tmpfile = files([impath filesep '*_peaks_*_structure..mat']);
    if ~isempty(tmpfile), delete(tmpfile{:}); end
end

%% Save to Excel %%
[day, tme] = bspm_timestamp;
outname = ['multi_table_I' num2str(intensity) '_S' num2str(separation) '_' day '_' tme '.xls'];
try 
    xlwrite(outname,combinedcell);
catch
    xlwrite(outname,combinedcell);
end

% cleanup
delete *_peaks_*structure.mat
delete *_peaks_*clusters.nii

















 
 
 
 
 
 
 
 
