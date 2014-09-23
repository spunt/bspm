function bspm_use_steffener(analysisdirs)
% BSPM_USE_STEFFENER
%
%   USAGE: [] = bspm_use_steffener(analysisdirs)
%
%   This is a personal modification of code by Dr. Jason Steffener which
%   combines beta images for canonical and derivate HRFs into a single
%   magnitude image (see below) for more description. 
%
%   ARGUMENTS
%       analysisdirs: analysis directory containing contrast images
%
% Text from Steffener Code:
% Inputs:
%           im1: full path to beta image corresponding to the primary
%              basis function
%           im2: full path to beta image corresponding to the secondary
%              basis function
%           Design: design matrix, e.g. SPM.xX.X which is the unfiltered
%               design matrix used stored in the SPM.mat file
%           ColumnsOfInterest: The two columns for which the contrast is
%               estimated for. The first column is expected to correspond to
%               the primary basis function.
%
% This program creates and executes the SPM command lines to calculate the magnitude
% image. It checks the design matrix to see if it is normalized.
%
% Written by: Jason Steffener
% js2746@columbia.edu
% date: November 3, 2009
% Bug Fixes:
% fix output file extensions

% ------------------------------- Copyright (C) 2014 -------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if ischar(analysisdirs), analysisdirs = cellstr(analysisdirs); end
nanalysis = length(analysisdirs);
for i = 1:nanalysis

    spmmat      = [analysisdirs{i} filesep 'SPM.mat'];
    tmp         = load(spmmat);
    nsess       = length(tmp.SPM.Sess); 
    xmatrix     = tmp.SPM.xX.X;
    betanam     = strcat(analysisdirs{i}, filesep, {tmp.SPM.Vbeta.fname}'); 
    xcon        = tmp.SPM.xCon;
    con_hdr     = {xcon.Vcon}';
    can_name    = regexprep({xcon.name}', ' - All Sessions', ''); 
    can_idx     = zeros(length(xcon), nsess);
    for c = 1:length(xcon)
        can_idx(c,:) = find(xcon(c).c); 
    end
    der_idx     = can_idx + 1;
    can_beta    = betanam(can_idx); 
    der_beta    = betanam(der_idx); 
    
    for c = 1:size(can_beta, 1)
        
        for s = 1:size(can_beta, 2)
        
            im1 = can_beta{c, s};
            im2 = der_beta{c, s};
            if size(can_beta, 2) > 1
                cname = ['MAG_' sprintf('%04d_SESS%d',c, s) '_' can_name{c} '.img'];
                cname = [analysisdirs{i} filesep cname];
                tmpn{s} = cname; 
            else
                cname = ['MAG_' sprintf('%04d',c) '_' can_name{c} '.img'];
                cname = [analysisdirs{i} filesep cname];
            end
            ccols = [can_idx(c,s) der_idx(c,s)];

            % BEGIN STEFFENER CODE
            Vin = spm_vol(im1);
            Vin(2) = spm_vol(im2);
            VOut = Vin(1);
            VOut.fname = cname;
            VOut.descrip = con_hdr{c}.descrip; 
            % Check Design matrix to determine if it is normalized
            tol = 0.0001;
            X = xmatrix(:,ccols);
            SSDesign = X'*X;
            %SSDesign = nX'*nX;
            if sum(diag(SSDesign) - [1; 1]) < tol
                % Design is normalized
            else
                nX(:,1) = X(:,1)./norm(X(:,1));
                nX(:,2) = X(:,2)./norm(X(:,2));
            end
            if SSDesign(2) > tol % Check off diagonal
                error('Design matrix is not orthogonalized');
            end
            % Create and execute the spm image calculation. 
            Weight1 = SSDesign(1,1);
            Weight2 = SSDesign(2,2);
            f = ['sqrt(i1.*i1.*' num2str(Weight1) ' + i2.*i2.*' num2str(Weight2) ')'];
            spm_imcalc(Vin,VOut,f);
            
        end
        
        if size(can_beta, 2) > 1
           
            h = spm_vol(char(tmpn));
            d = spm_read_vols(h);
            hout = h(1); 
            hout.fname = regexprep(hout.fname, 'SESS1_', '');
            spm_write_vol(hout, mean(d, 4)); 
            
        end

    end

    
end






 
 
 
 
