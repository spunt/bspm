function bspm_auto_check_reg(imcell,subname,xpos,imname)
% BSPM_AUTO_CHECK_REG Identifies scans and saves a nuisance regressor file
% 
% USAGE: bspm_auto_check_reg(imcell,subname,xpos,imname)
%
%  ARGUMENTS
%   imcell = nsubs x nimages cellarray 
%  

% ---------------------------- Copyright (C) 2014 ----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, error('USAGE: bspm_auto_check_reg(imcell,subname,xpos,imname)'); end;
if nargin<3, xpos = 0; end
if nargin<4, 
    imname = cell(size(imcell,2),1); 
    for i = 1:length(imname), 
        [path, name, ext] = fileparts(imcell{1,i});
        [path, imname{i}, ext] = fileparts(path);
    end 
end
for i = 1:size(imcell,1)
    cim = imcell(i,:);
    csub = subname{i};
    for p = 1:length(xpos)
        bspm_checkreg(cim, imname, [xpos(p) 40 0], [csub ',a =  x=' num2str(xpos(p))]);
        saveas(gcf, sprintf('%s_x=%d.pdf', csub, xpos(p)), 'pdf');
    end
    close all;
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
