function peaks = bspm_get_peaks(image,intensity,cluster,separation,notstat)
% BSPM_SAVE_TABLE
%
%   USAGE: bspm_get_peaks(image,intensity,cluster,separation)
%
%   ARGUMENTS
%       image:      filename of input statistic image
%       intensity:  intensity threshold to use
%       cluster:    cluster extent threshold to use
%       separation: minimum peak separation (in mm)
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

if nargin<5, notstat = 0; end
if nargin<4, error('USAGE: bspm_get_peaks(image,intensity,cluster,separation)'); end

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
peaks1 = voxels{1};
peaks1(:,6:7) = [];

peaks.labels = labels;
peaks.coords = peaks1(:,3:5);
peaks.k = peaks1(:,1);
peaks.tstat = peaks1(:,2);


% cleanup
delete *_peaks_*structure.mat
delete *_peaks_*clusters.nii

 
 
 
 
