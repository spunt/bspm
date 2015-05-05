function cmd = bfsl_ica_aroma(infile, TR, outdir)
% BFSL_ICA_AROMA Run ICA-AROMA using FSL
%
%  USAGE: bfsl_ica_aroma(infile, TR, outdir)
% __________________________________________________________________________
%  INPUTS
%	infile:     4D timeseries to denoise
%	TR:         in seconds [optional]
%	outdir:     output directory [default = same as in]
%
%
% usage: ICA_AROMA.py [-h] -o OUTDIR [-i INFILE] [-mc MC] [-a AFFMAT] [-w WARP]
%                     [-m MASK] [-f INFEAT] [-tr TR] [-den DENTYPE] [-md MELDIR]
%                     [-dim DIM]
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-24
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 1, disp('USAGE: bfsl_ica_aroma(infile, TR, outdir)'); return; end
if iscell(infile), infile = char(infile); end
if nargin < 2, TR = []; end
if nargin < 3, outdir = fileparts(infile); end

% | - Configure Path
aromadir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'thirdparty', 'ICA-AROMA');  
icaaroma = fullfile(aromadir, 'ICA_AROMA.py'); 
rpfile      = char(files(fullfile(fileparts(infile), 'rp*txt')));
if isempty(TR)
    cmd = sprintf('python2.7 %s -i %s -o %s -mc %s &', icaaroma, infile, outdir, rpfile);
else
    cmd = sprintf('python2.7 %s -i %s -o %s -mc %s -tr %d &', icaaroma, infile, outdir, rpfile, TR);
end
if nargin==0
    system(cmd);
end