function bspm_use_steffener_old(analysisdirs)
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

    spmmat = [analysisdirs{i} filesep 'SPM.mat'];
    tmp = load(spmmat);
    xmatrix = tmp.SPM.xX.X;
    connam = tmp.SPM.xX.name;
    betanam = {tmp.SPM.Vbeta.fname};
    xcon    = tmp.SPM.xCon;
    con_hdr = {xcon.Vcon}';
    for r = 1:length(tmp.SPM.Sess)
        string = sprintf('Sn\\(%d\\) ', r);
        connam = regexprep(connam,string,'');
    end
    can_idx = cellstrfind(connam,'bf(1)');    %  canonical 
    der_idx = cellstrfind(connam,'bf(2)');    % derivative
    can_beta = betanam(can_idx);
    can_name = connam(can_idx);
    for c = 1:length(can_name)
        n = can_name{c};
        n(strfind(n,'*'):end) = [];
        can_name{c} = n;
    end
    der_beta = betanam(der_idx);
    der_name = connam(der_idx);
    
    for c = 1:length(can_beta)
        
        im1 = [analysisdirs{i} filesep can_beta{c}];
        im2 = [analysisdirs{i} filesep der_beta{c}];
        cname = ['mag_' sprintf('%04d',c) '_' can_name{c} '.img'];
        cname = [analysisdirs{i} filesep lower(cname)];
        ccols = [can_idx(c) der_idx(c)];
        
        % BEGIN STEFFENER CODE
        Vin = spm_vol(im1);
        Vin(2) = spm_vol(im2);
        VOut = Vin(1);
        VOut.fname = cname;
        VOut.descrip = [con_hdr{c}.descrip ' Magnitude Image']; 
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

    
end






 
 
 
 
