function bspm_con2tstat(con)
% BSPM_CON2TSTAT
%
%   USAGE: [tstat df] = bspm_con2tstat(con)
%       
%       con  =  array of images OR wildcard pattern for finding them
%       nowrite = flag to not write spmT image (default = 0)
%

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1, mfile_showhelp; return; end
if ischar(con), con = cellstr(con); end

for i = 1:length(con)
    
    %% find SPM.mat and ResMS.img
    [p n e] = fileparts(con{i});
    try
        load([p filesep 'SPM.mat']);
    catch
        error('No SPM.mat in %s', p); 
    end
    try
        tmpresms = [p filesep 'ResMS.img'];
        if ~exist(tmpresms, 'file')
            tmpresms = [p filesep 'ResMS.nii'];
        end
        h = spm_vol(tmpresms);
        RESMS = spm_read_vols(h);
    catch
        error('No ResMS.img found in %s', p); 
    end
    DF = SPM.xX.erdf;
    BCOV = SPM.xX.Bcov;
    hdr = spm_vol(con{i});
    DATA = spm_read_vols(hdr);
    descrip = regexprep(hdr.descrip,'SPM contrast - ','');
    tmpidx = 1; 
    connum = 1; 
    conname = strtrim(descrip(tmpidx(1)+1:end));
    conname = regexprep(conname, ' - All Sessions','');
    CONWEIGHT = SPM.xCon(connum).c;
    VcB = CONWEIGHT'*BCOV*CONWEIGHT;
    T = DATA./sqrt(RESMS*VcB);
    thdr = hdr;
    thdr.descrip = sprintf('SPM{T_[%2.1f]} - %s', DF, conname);
    thdr.fname = sprintf('%s/spmT_000%d_%s.nii', p, connum, conname);
    spm_write_vol(thdr, T);

end

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
