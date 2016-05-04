function outname = bspm_coords2roi(xyz,roi,labels)
% BSPM_COORDS2ROI
%
%   USAGE: bspm_coords2roi(ref,xyz,roi)
%
%   ARGUMENTS
%       xyz = coordinates (rows are ROIs)
%       roi.shape = 'Sphere' or 'Box'
%       roi.size = radius of ROI
%       
% Uses code authored by:
% Dr. Robert Welsh (SimpleROIBuilder.m)
% Drs. Donald McLaren & Aaron Schultz (peak_nii.m)

% -------------- Copyright (C) 2014 --------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_12_2014

if nargin<2, mfile_showhelp; return; end
if nargin<3, labels = repmat({''},size(xyz,1),1); end
spmdir = fileparts(which('spm'));
ref = fullfile(spmdir, 'templates', 'EPI.nii'); 
if ischar(labels), labels = cellstr(labels); end
refhdr = spm_vol(ref);
roihdr = refhdr;
roihdr.pinfo = [1;0;0];
roipath = pwd;
[R,C,P]  = ndgrid(1:refhdr.dim(1),1:refhdr.dim(2),1:refhdr.dim(3));
RCP      = [R(:)';C(:)';P(:)'];
clear R C P
RCP(4,:) = 1;
XYZmm    = refhdr.mat(1:3,:)*RCP;   
Q          = ones(1,size(XYZmm,2));
nroi = size(xyz,1);
for i = 1:nroi
    
    mm = xyz(i,:)';
    cROI = zeros(roihdr.dim);
    cHDR = roihdr;
    roidescrip = ['ROI_' roi.shape num2str(roi.size) '_' num2str(round(xyz(i,1))) '_' num2str(round(xyz(i,2))) '_' num2str(round(xyz(i,3))) '_' labels{i}];
    cHDR.fname = [roipath filesep roidescrip '.nii'];
    cHDR.descrip = roidescrip;
    switch roi.shape
        case 'Sphere'
        j = find(sum((XYZmm - mm*Q).^2) <= roi.size^2);
        case 'Box'
        j      = find(all(abs(XYZmm - mm*Q) <= [roi.size roi.size roi.size]'*Q/2));
    end
    cROI(j) = 1;
    outname = cHDR.fname; 
    spm_write_vol(cHDR,cROI);
    fprintf('ROI file created: %s\n', [roidescrip '.nii']);
    
end
    


 
 
 
 
