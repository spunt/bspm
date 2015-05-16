function [] = bspm_imcalc_old(in, outname, operation)
% BSPM_IMCALC  Calculate Image based on Input Image(s)
%
%   USAGE: bspm_imcalc(in, outname, operation)
%       
%       in  =  array of images to smooth (full path)
%       outname = name for output image
%       operation = string specifying operation to apply
%
%           LOGICAL OPERATORS
%               >, <, >=, <=, ==
%               
%           NON-LOGICAL OPERATORS (ACROSS IMAGES)
%               'sum'   - sum across images
%               'prod'  - produce across images
%               'mean'  - mean across images
%               'median' - median across images
%               'std'   - std across images
%               'var'   - var across images
%               'min'   - min across images
%               'max'   - max across images
%               'diff'  - contrast across images (2 only)
%
%           NON-LOGICAL OPERATORS (SINGLE IMAGE)
%               'sqrt'          - square root of image
%               'negative'      - negative of image
%               'nan2zero'      - convert NaNs to 0s
%               'zero2nan'      - convert 0s to NaNs
%               'zscore'        - convert t-image to z-image
%               'prctile'       - convert to percentile
%               'prctilesym'    - convert to percentile for +/- separately
%               'fill'          - uses IMFILL to fill holes in volume
%
%           SPECIAL OPERATORS
%               'colorcode' - combine and colorcode multiple images
%
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2013.03.20
% CREDIT: Loosely based on functionality of spm_imcalc_ui.m (SPM8)

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, disp('USAGE: bspm_imcalc(in, outname, operation)'); return, end

% | check variable formats
if ischar(in), in = cellstr(in); end
if iscell(operation), operation = char(operation); end
    
% | read in data
[im, hdr] = bspm_read_vol(in);
if length(hdr)>1, hdr = hdr(1); end

% | apply the operation
tag = 0;
switch lower(operation)
    case {'prod'}
        outim = eval([operation '(im,4);']);
    case {'sum', 'mean', 'median'}
        operation = ['nan' operation];
        outim = eval([operation '(im,4);']);
    case {'min', 'max', 'std', 'var'}
        operation = ['nan' operation];
        outim = eval([operation '(im,[],4);']);
    case {'negative'}
        outim = -1*im;
    case {'sqrt'}
        outim = sqrt(im);
    case {'diff'}
        outim = -1*(eval([operation '(im,1,4);']));
    case {'nan2zero'}
        outim = im; outim(isnan(outim)) = 0; tag = 1;
    case {'zero2nan'}
        outim = im; outim(outim==0) = NaN; tag = 1;
    case {'colorcode'}
        outim = im(:,:,:,1);
        for i = 1:length(in)
            tmpim = im(:,:,:,i);
            idx = tmpim(:)>0;
            outim(idx) = i;
        end
    case {'zscore'}
        i1 = strfind(hdr.descrip,'[');
        i2 = strfind(hdr.descrip,']');
        df = str2num(hdr.descrip(i1+1:i2-1));
        outim = im; 
        outim(abs(outim) > 0) = bspm_t2z(outim(abs(outim) > 0),df);
        outim(outim==Inf) = max(outim(outim~=Inf))*1.01;
    case {'prctile'}
        outim = im; 
        outim(abs(outim) > 0) = 100*(tiedrank(outim(abs(outim) > 0)) ./ sum(abs(outim(:)) > 0));
    case {'prctilesym'}
        outim = im;
        posidx = find(outim > 0); 
        negidx = find(outim < 0); 
        pos = 100*(tiedrank(outim(posidx)) ./ length(posidx));
        neg = -100*(tiedrank(abs(outim(negidx))) ./ length(negidx));
        outim(posidx) = pos;
        outim(negidx) = neg;
    case {'imfill'}
        outim = imfill(im,6,'holes');
    otherwise
        outim = eval(['im' operation]);
end

% construct outname and write image
if ~tag
    tmp = hdr.descrip;
    idx = regexp(tmp,'SPM{T','ONCE');
    if ~isempty(idx)
        tmp(regexp(tmp,'-')+2:end) = [];
        hdr.descrip = [tmp upper(operation) ' IMAGE'];
    else
        hdr.descrip = [upper(operation) ' IMAGE'];
    end
end
hdr.fname = outname; 
spm_write_vol(hdr, outim); 
% bnii_write(outim, hdr, 'outname', outname)    
 
 
 
 
 
 
 
 
