function bspm_append_descrip(in, deletepat)
% BSPM_APPEND_DESCRIP
%
%   USAGE: bspm_append_descrip(in, deletepat)
%       
%       in  =  array of spmT* or con* images to rename
%       deletepat = parts of descrip to delete
%

% ------------------ Copyright (C) 2014 ------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
defaultpat = {'spm_spm:beta \(\d\d\d\d\)\s-\sSn\(\d\)'
'\*.*$'
'^SPM\{T_\[\d+\.\d\]\}\s-\scontrast\s\d+:'
'^Contrast\s\d+:'
'-\sAll\sSessions'};
if nargin<1, mfile_showhelp; return; end
if nargin<2, deletepat = []; end
if ischar(in), in = cellstr(in); end
if ischar(deletepat), deletepat = cellstr(deletepat); end
if size(deletepat, 2) > 1, deletepat = deletepat'; end
deletepat = [defaultpat; deletepat];
[pth, fname, fext] = cellfun(@fileparts, in, 'Unif', false); 
nimg = length(in);
hdr = spm_vol(char(in));
des = {hdr.descrip}';
outname = des; 
for i = 1:length(deletepat), outname = strtrim(regexprep(outname, deletepat{i}, '')); end
out = fullfile(pth, strcat(fname, '_', outname, fext));
cellfun(@copyfile, in, out); 

end





% for i = 1:nimg
%     
% 
%     % get name
%     tmp = hdr.descrip;
%     p1 = strfind(tmp,':') + 1;
%     tmpidx = strfind(tmp,'All Sessions');
%     if ~isempty(tmpidx)
%         connam = tmp(p1:tmpidx-3);
%     elseif strfind(tmp,'beta')
%         tmpidx = strfind(tmp,' - ');
%         tmpidx2 = strfind(tmp,'*');
%         connam = tmp(tmpidx(end)+9:tmpidx2(end)-1);
%     else
%         connam = tmp(p1:end);
%     end
%     connam = strtrim(connam);
%     % write new image
%     oldname = hdr.fname;
%     [p n e] = fileparts(oldname);
%     hdr.fname = [p filesep n '_' connam '.' format];
%     if rmTAG
%         delete([n '*']);
%     end
%     spm_write_vol(hdr, img);
%     
% end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
