function matlabbatch = bspm_3dto4d(in, outname)
% BSPM_3DTO4D
%
% USAGE: matlabbatch = bspm_3dto4d(in, outname)
%
% ARGUMENTS
%   in: array of 3D volumes to convert (cell or char)
%   out: name to give 4D volume [default = 4D.nii]
%
% Bob Spunt, November 18, 2012
% 

% ---------------- Copyright (C) 2014 ----------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin < 1, mfile_showhelp; return; end
if ischar(in), in = cellstr(in);end
if nargin < 2
    [p,n,e] = cellfun(@fileparts, in, 'unif', false);
    n       = regexp(n, '\W', 'split');
    nelem   = length(n{1});  
    n       = [n{:}]; 
    n       = reshape(n, nelem, length(p))';
    for i = 1:nelem, nu(i) = length(unique(n(:,i))); end
    n       = n(1, nu==1);
    n       = strcat(n, '-');
    n       = strcat(n{:}, '4D.nii');
    outname = fullfile(p{1}, n); 
end
matlabbatch{1}.spm.util.cat.vols    = in; 
matlabbatch{1}.spm.util.cat.name    = outname;
matlabbatch{1}.spm.util.cat.dtype   = 0;

% | Run job (only if no output arguments requested)
% | =======================================================================
if nargout==0, spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end

% 
% V    = spm_vol(char(in));
% ind  = cat(1,V.n);
% N    = cat(1,V.private);
% mx   = -Inf;
% mn   = Inf;
% for i=1:numel(V),
%     dat      = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
%     dat      = dat(isfinite(dat));
%     mx       = max(mx,max(dat(:)));
%     mn       = min(mn,min(dat(:)));
% end;
% dat(dat==0) = NaN;
% sf         = max(mx,-mn)/32767;
% ni         = nifti;
% ni.dat     = file_array(outname,[V(1).dim numel(V)],'INT16-BE',0,sf,0);
% ni.mat     = N(1).mat;
% ni.mat0    = N(1).mat;
% ni.descrip = '4D image';
% create(ni);
% for i=1:size(ni.dat,4),
%     ni.dat(:,:,:,i) = N(i).dat(:,:,:,ind(i,1),ind(i,2));
%     spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
% end;
%  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
