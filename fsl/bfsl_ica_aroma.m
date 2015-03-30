function bfsl_ica_aroma(in,TR,outdir)
% BFSL_ICA_AROMA Run ICA-AROMA using FSL
%
%  USAGE: bfsl_ica_aroma(in,TR,outdir)
% __________________________________________________________________________
%  INPUTS
%	infile:     4D timeseries to denoise
%	TR:         in seconds
%	outdir:     output directory
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
if nargin < 1, disp('USAGE: bfsl_ica_aroma(in, TR, outdir)'); return; end
if iscell(in), in = char(in); end
if nargin < 2, TR = []; end
if nargin < 3, outdir = fileparts(in); end
rpfile = char(files(fullfile(fileparts(in), 'rp*txt'))); 
icaaroma = '/Users/bobspunt/Desktop/Dropbox/Bob/Matlab/_functions_/git/others/ICA-AROMA/ICA_AROMA.py';
if isempty(TR)
    cmd = sprintf('python2.7 %s -i %s -o %s -mc %s', icaaroma, in, outdir, rpfile);
else
    cmd = sprintf('python2.7 %s -i %s -o %s -mc %s -tr %d', icaaroma, in, outdir, rpfile, TR);
end
system(cmd);