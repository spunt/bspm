function bspm_save_table(image,intensity,cluster,separation,doinv,notstat)
% BSPM_SAVE_TABLE
%
%   USAGE: bspm_save_table(image,intensity,cluster,separation,doinv,notstat)
%
%   ARGUMENTS
%       image:      filename of input statistic image
%       intensity:  intensity threshold to use
%       cluster:    cluster extent threshold to use
%       separation: minimum peak separation (in mm)
%       doinv:      flag to do inverse (default = 1)
%       notstat:    flag for image that isn't statistical
%      
% Created April 8, 2013 - Bob Spunt
% Uses code authored by:
% Drs. Donald McLaren & Aaron Schultz (peak_nii.m) 

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, separation = 20; end
if nargin<4, mfile_showhelp; return; end
if nargin<5, doinv = 1; end
if nargin<6, notstat = 0; end

basedir = pwd; 

% make sure image name is character array
if iscell(image), image = char(image); end

% default structure input to peak_nii
peaknii.thresh = intensity;
peaknii.cluster = cluster;
peaknii.out = '';
peaknii.sign = 'pos';
if notstat
    peaknii.type = 'none';
else
    peaknii.type = 'T';
end
peaknii.voxlimit = [];
peaknii.separation = separation;
peaknii.SPM = 1;
peaknii.conn = [];
peaknii.mask = [];
peaknii.df1 = [];
peaknii.df2 = [];
peaknii.nearest = 1;
peaknii.label = 'HarvardOxford_cortex';

% run peak_nii
[voxels] = peak_nii(image,peaknii);

% clean up data
labels = voxels{2};
labels = regexprep(labels,',','');
peaks = voxels{1};
peaks(:,6:7) = [];

% make inverse if appropriate
[impath imname e] = fileparts(image);
if doinv
    peaknii.sign = 'neg';
    [voxels] = peak_nii(image,peaknii);
    invlabels = voxels{2};
    invlabels = regexprep(invlabels,',','');
    invpeaks = voxels{1};
    invpeaks(:,6:7) = [];
    invpeaks(:,2) = abs(invpeaks(:,2));
end

% write table
if doinv
    headers1 = {'' '' '' '' '' 'MNI Coordinates' '' ''};
    headers2 = {'Contrast' '' '' '' '' '' '' ''};
    headers3 = {'' 'Region Name' 'L/R' 'Extent' 't-value' 'x' 'y' 'z'};
    headers4 = {'Positive' '' '' '' '' '' '' ''};
    headers5 = {'Negative' '' '' '' '' '' '' ''};
    [p imfolder e] = fileparts(impath);
    outname = ['save_table_' imfolder '_' imname '_I' num2str(intensity) '_C' num2str(cluster) '_S' num2str(separation) '.xls'];
    side = labels;
    side(:) = {'L'};
    side(peaks(:,3) > 0) = {'R'};
    side(peaks(:,3)==0) = {'L/R'};
    pad = cell(size(side));
    datacell = [pad labels side num2cell(peaks)];
    side = invlabels;
    side(:) = {'L'};
    side(invpeaks(:,3) > 0) = {'R'};
    side(invpeaks(:,3)==0) = {'L/R'};
    pad = cell(size(side));
    invdatacell = [pad invlabels side num2cell(invpeaks)];
    rowpad = cell(size(headers1));
    allcell = [rowpad; headers1; headers2; headers3; rowpad; headers4; datacell; rowpad; headers5; invdatacell];
    xlwrite(outname,allcell);
else
    headers1 = {'' '' '' 'MNI Coordinates' '' ''};
    headers2 = {'Region Name' 'Extent' 't-value' 'x' 'y' 'z'};
    [p imfolder e] = fileparts(impath);
    outname = ['save_table_' imfolder '_' imname '_I' num2str(intensity) '_C' num2str(cluster) '_S' num2str(separation) '.xls'];
    datacell = [labels num2cell(peaks)];
    allcell = [headers1; headers2; datacell];
    xlwrite(outname,allcell);
end

% cleanup
cd(impath)
delete *_peaks_*structure.mat
delete *_peaks_*clusters.nii
cd(basedir)

 
 
 
 
