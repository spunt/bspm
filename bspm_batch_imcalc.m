function [] = bspm_batch_imcalc(in, prefix, operation)
% BSPM_BATCH_IMCALC  Simple wrapper for bspm_imcalc
%
%   USAGE: bspm_batch_imcalc(in, prefix, operation)
%       
%       in  =  cell array of images (full path)
%       prefix = prefix for each output image
%       operation = string specifying operation to apply
%          
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2013.03.20
% CREDIT: Loosely based on functionality of spm_imcalc_ui.m (SPM8)

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, disp('USAGE: bspm_batch_imcalc(in, prefix, operation)'); return, end

% check variable formats
if ischar(in), in = cellstr(in); end
if iscell(operation), operation = char(operation); end
    
% read in data
for i = 1:length(in)
    
    cim = in{i};
    name = [];
    if ~isempty(prefix)
        [p,tmpname,ext] = fileparts(in{i});
        name = [name '_' tmpname];
        outname = [path filesep prefix name ext];
    else
        outname = cim;
    end
    bspm_imcalc(in{i},outname,operation);

end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
