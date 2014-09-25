function bspm_add_connam(in, rmTAG, format)
% BSPM_ADD_CONNAM
%
%   USAGE: bspm_add_connam(in, rmTAG)
%       
%       in  =  array of spmT* or con* images to rename
%       rmTAG = option to delete original (default = 0)
%       format = of output image, 'img' or 'nii' (default)
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<3, format = 'nii'; end
if nargin<2, rmTAG = 0; end
if nargin<1, disp('USAGE: bspm_add_connam(in, rmTAG)'); return; end
if ~iscell(in) & strfind(in,'*'); in = files(in); end

% make sure image names are cell arrays of strings
if strfind(format,'.'), format = regexprep(format, '\.', ''); end
if ischar(in), in = cellstr(in); end

% loop over volumes, add name, then write new volume
for i = 1:length(in)
    
    hdr = spm_vol(in{i});
    img = spm_read_vols(hdr);
    % get name
    tmp = hdr.descrip;
    p1 = strfind(tmp,':') + 1;
    tmpidx = strfind(tmp,'All Sessions');
    if ~isempty(tmpidx)
        connam = tmp(p1:tmpidx-3);
    elseif strfind(tmp,'beta')
        tmpidx = strfind(tmp,' - ');
        tmpidx2 = strfind(tmp,'*');
        connam = tmp(tmpidx(end)+9:tmpidx2(end)-1);
    else
        connam = tmp(p1:end);
    end
    connam = strtrim(connam);
    % write new image
    oldname = hdr.fname;
    [p n e] = fileparts(oldname);
    hdr.fname = [p filesep n '_' connam '.' format];
    if rmTAG
        delete([n '*']);
    end
    spm_write_vol(hdr, img);
    
end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
