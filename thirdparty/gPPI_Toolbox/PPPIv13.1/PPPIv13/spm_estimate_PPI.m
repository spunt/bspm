function SPM=spm_estimate_PPI(Subject,SPM,Regions,Method,Analysis,Contrast,ConcatR,outdir,preservevarcorr,wb)
%Estimates 1st level PPI model
%   Subject is the subject number 
%   SPM is the 1st level activity model, used to pull values from
%   Regions are the PPI seed regions
%   Method is 'cond' or 'trad'
%   Analysis is either 'psy', 'phy', or 'psyphy'
%   Contrast is a computational method variable:
%       0 not to estimate any contrasts;
%       1 to estimate contrasts;  
%       2 to only use PPI txt file for 1st level (not recommended); 
%       3 to only use PPI txt file for 1st level and estimate contrasts (not recommended);
%           NOTE: 2&3 are not recommended as they potentially do not include
%                 all tasks effects in the mode. Use at your own risk.
%           NOTE: Not required, defaults to act like 0.
%
% License:
%   Copyright (c) 2011,2012 Donald G. McLaren and Aaron Schultz
%   All rights reserved.
%
%    Redistribution, with or without modification, is permitted provided that the following conditions are met:
%    1. Redistributions must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%    2. All advertising materials mentioning features or use of this software must display the following acknowledgement:
%        This product includes software developed by the Harvard Aging Brain Project.
%    3. Neither the Harvard Aging Brain Project nor the
%        names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%    4. You are not permitted under this Licence to use these files
%        commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all 	
%        or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use 	
%        of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third 	
%        party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale 
%        or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received.
%
%   THIS SOFTWARE IS PROVIDED BY DONALD G. MCLAREN (mclaren@nmr.mgh.harvard.edu) AND AARON SCHULTZ (aschultz@nmr.mgh.harvard.edu)
%   ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
%   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Last modified on 3/14/2012 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%   Medical School
%
%   Update on 3/14/2012
%   Added the option to use temporal/dispersion derivatives
%   Update on 11/28/2012
%   Added option for specifying output directory
%
%   Update in March 2013
%   Added option for using pre-loaded data to increase speed for multiple ROI.
%   Added option to use existing V structure from first-level model.

%   Update in April 2014
%   Added option FSFAST and rWLS
%   Added the ability to switch between spm_spm and spm_spm_WB internally

%% Check parameters and setup variables
if exist(SPM,'file')
    load(SPM)
elseif exist(SPM,'var')
else
    disp('Program will now exit. No SPM.mat was detected.')
    return;
end
% Directory
[region1,region2]=strtok(Regions);
region2=region2(2:end);
if strcmp(Analysis,'psyphy')
    A='PPPI';
else
    A='PPI';
end
oSPMdir=SPM.swd; % Used for storing datamatrix for repeated ROIs of the same subject
if isempty(region2)
   if exist('outdir','var') && ~isempty(outdir)
       SPMPPI.swd=[outdir filesep A '_'  region1 filesep];
   else
       SPMPPI.swd=[SPM.swd filesep A '_'  region1 filesep];
   end
else
    if exist('outdir','var') && ~isempty(outdir)
        SPMPPI.swd=[outdir filesep A '_' region1 '_' region2 filesep];
    else
        SPMPPI.swd=[SPM.swd filesep A '_' region1 '_' region2 filesep];
    end
end

if exist([SPMPPI.swd 'SPM.mat'],'file')==2
    SPMPPI2=load([SPMPPI.swd 'SPM.mat']);
    if isfield(SPMPPI2.SPM,'VM') && (exist(SPMPPI2.SPM.VM.fname,'file')==2 || exist([SPMPPI.swd SPMPPI2.SPM.VM.fname],'file')==2)
        estimate=0;
        SPMPPI=SPMPPI2; clear SPMPPI2
    else
        estimate=1;
    end
else
    estimate=1;
end

if estimate==1
    % PPI File Extension
    try
        if Contrast==2 || Contrast==3
            ext='.txt';
        else
            ext='.mat';
        end
    catch
        ext='.mat';
    end

    % Get needed fields from task model
    SPMPPI.xBF=SPM.xBF;
%     deriv=0; % do not include derivative terms from first level task model
%     %deriv=1; % include derivative terms from first level task model
%     if deriv==0;
%         SPMPPI.xBF.name='hrf';
%         SPMPPI.xBF.order=1;
%         SPMPPI.xBF.bf=SPMPPI.xBF.bf(:,1:SPMPPI.xBF.order);
%     end
    SPMPPI.xY=SPM.xY;
    SPMPPI.nscan=SPM.nscan;
    SPMPPI.SPMid=SPM.SPMid;
    SPMPPI.xVi.form=SPM.xVi.form;
    if ischar(SPM.xVi.form) && strcmp(SPM.xVi.form,'i.i.d')
        SPMPPI.xVi.form='none';
    end 
    SPMPPI.xGX.iGXcalc=SPM.xGX.iGXcalc;
    SPMPPI.xGX.sGXcalc=SPM.xGX.sGXcalc;
    SPMPPI.xGX.sGMsca=SPM.xGX.sGMsca;
    for i=1:numel(SPMPPI.nscan)
        SPMPPI.xX.K(i).HParam = SPM.xX.K(i).HParam;
        SPMPPI.xX.K(i).RT = SPM.xX.K(i).RT;
    end

    for z=1:numel(SPM.Sess)
        if strcmp(ext,'.txt')
            try
                load([Subject '_' region1 '_and_' region2 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            catch
                load([Subject '_' region1 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            end
            SPMPPI.Sess(z).C.C=[OUT.P.C OUT.PPI.C OUT.Y.C OUT.C.C];
            SPMPPI.Sess(z).C.name=[OUT.P.name OUT.PPI.name OUT.Y.name OUT.C.name];
            return;
        else
            try
                load([Subject '_' region1 '_and_' region2 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            catch
                load([Subject '_' region1 '_session' num2str(z) '_' Method '_' A '_regressors' ext]);
            end
            if strcmpi(Method,'trad')
                SPMPPI.Sess(z).U=[];
                SPMPPI.Sess(z).C.C=[OUT.P.C OUT.PPI.C OUT.Y.C OUT.C.C];
                SPMPPI.Sess(z).C.name=[OUT.P.name OUT.PPI.name OUT.Y.name OUT.C.name];
            else
                SPMPPI.Sess(z).U=SPM.Sess(z).U;
                SPMPPI.Sess(z).C.C=[OUT.PPI.C OUT.Y.C OUT.C.C];
                SPMPPI.Sess(z).C.name=[OUT.PPI.name OUT.Y.name OUT.C.name];
            end
        end
    end
    
    % Concatenated PPI, concatenates the runs, but preserves the filtering
    % and uses AR1 filtering from first level (see below). Beta Version.
    if exist('ConcatR','var') && ConcatR==1 
        if strcmpi(Method,'trad')
            error('ConcatR option not valid with traditional PPI')
        else
            % Make K compatible with design, will revert back to per
            % session filtering after design creation
            SPMPPI.xX.K(1).row=SPM.xX.K(1).row;
            for ii=2:numel(SPMPPI.nscan)
                SPMPPI.xX.K(1).row=[SPMPPI.xX.K(1).row SPM.xX.K(ii).row];   
            end
            SPMPPI.nscan=sum(SPMPPI.nscan);
            SPMPPI.xX.K(2:end)=[]; %remove K for subsequent runs, since this is a single run model
            Sess.U=[];
            
            %Task Regressors
            tasks={};
            for ii=1:numel(SPMPPI.Sess)
                for jj=1:numel(SPMPPI.Sess(ii).U)
                    tasks=[tasks SPMPPI.Sess(ii).U(jj).name];
                end
            end
            [Tasks]=unique(tasks);
            clear tasks
            for tt=1:numel(Tasks)
                ind=contains([' ?' Tasks{tt}],SPM.xX.name);
                Sess.C.C(:,tt)=sum(SPM.xX.X(:,ind),2);
                Sess.C.name{tt}=Tasks{tt};
            end
            cols=tt;
            
            %PPI Regressors
            mast_ind=cell(numel(SPMPPI.Sess),1);
            for tt=1:numel(Tasks)
                Sess.C.C(:,cols+tt)=0;
                for rr=1:numel(SPMPPI.Sess)
                    ind=contains(['PPI_' Tasks{tt}],SPMPPI.Sess(rr).C.name);
                    mast_ind{rr}=[mast_ind{rr} ind];
                    if ~isempty(ind)
                       Sess.C.C(SPM.Sess(rr).row,cols+tt)=SPMPPI.Sess(rr).C.C(:,ind);
                    end
                end
                Sess.C.name{cols+tt}=['PPI_' Tasks{tt}];
            end
            
            % Seed Regressor
            ind=contains('seedtc',SPMPPI.Sess(1).C.name);
            seed{1}=SPMPPI.Sess(1).C.name{ind(1)};
            if strcmpi(Analysis,'phy') || strcmpi(Analysis,'psyphy')
                seed{2}=SPMPPI.Sess(1).C.name{ind(2)};
            end
            for tt=1:numel(seed)
                Sess.C.C(:,end+1)=0;
                cols=size(Sess.C.C,2);
                for rr=1:numel(SPMPPI.Sess)
                    ind=contains(seed{tt},SPMPPI.Sess(rr).C.name);
                    mast_ind{rr}=[mast_ind{rr} ind];
                    Sess.C.C(SPM.Sess(rr).row,cols)=SPMPPI.Sess(rr).C.C(:,ind);
                end
                Sess.C.name{end+1}=seed{tt};
            end
            
            % Other Regressors
            for rr=1:numel(SPMPPI.Sess)
                ind=setdiff(1:numel(SPMPPI.Sess(rr).C.name),mast_ind{rr});
                Sess.C.C(:,end+1:end+numel(ind))=0;
                Sess.C.C(SPM.Sess(rr).row,end-numel(ind)+1:end)=SPMPPI.Sess(rr).C.C(:,ind)-repmat(mean(SPMPPI.Sess(rr).C.C(:,ind)),size(SPMPPI.Sess(rr).C.C,1),1);
                Sess.C.name=[Sess.C.name SPMPPI.Sess(rr).C.name{ind}];
            end
            
            %Constants (for runs)
            for rr=2:numel(SPMPPI.Sess)
                Sess.C.C(:,end+1)=0;
                Sess.C.C(SPM.Sess(rr).row,end)=1;
                Sess.C.name{end+1}=['Constant_' num2str(rr)];
            end
            SPMPPI.SessC=SPM.Sess;
            SPMPPI.Sess=Sess;
        end
    end
        
    %% Estimate PPI 1st Level Model
    try
        cd(SPMPPI.swd)
    catch
        mkdir(SPMPPI.swd)
        cd(SPMPPI.swd)
    end
    save SPM SPMPPI

    % Delete any existing files
    delete beta_00*
    delete ResMS.*
    delete RPV.*
    delete mask.*

    % Make design and estimate, rename to SPM
    SPM1=SPM; clear SPM
    SPM=SPMPPI; clear SPMPPI
    SPM = spm_fmri_spm_ui(SPM);
    
    disp('estimate_PPI.m')
    SPM.xM.T(:) = -Inf;  %% disable threshold masking
    SPM.xM.TH(:)=SPM1.xM.TH(:);
    try
        V = spm_vol(SPM1.VM.fname);
    catch
        try
            V = spm_vol([SPM1.swd filesep SPM1.VM.fname]);
        catch
            disp(['Mask file cannot be found: ' SPM1.VM.fname])
            error('Mask file cannot be found.')
        end
    end
    SPM.xM.VM = V;
    % Preserve HP filtering; Use AR1 filtering from First Level Model - DO
    % NOT RECOMPUTE VARIANCE CORRECTION.
    if exist('ConcatR','var') && ConcatR==1
        if strcmpi(Method,'trad')
            error('ConcatR option not valid with traditional PPI')
        else
            SPM.xX.K=SPM1.xX.K;
            SPM.xVi=SPM1.xVi;
        end
    end
    if exist('preservevarcorr','var') && preservevarcorr==1
        SPM.xVi=SPM1.xVi;
    end
    if exist('preservevarcorr','var') && preservevarcorr==2
        SPM.xVi=SPM1.xVi;
        SPM.xX.W=SPM1.xX.W;
    end
    if ~isempty(strfind(SPM.swd,['.lh' filesep])) || ~isempty(strfind(SPM.swd,['.rh' filesep]))
        SPM.xVol.FWHM=[];
        SPM.xVol.VRpv=[];
        SPM.xVol.R=[];
    end
    if exist('wb','var') && wb==1
        SPM=spm_spm_WB(SPM,[oSPMdir filesep 'datamatrixforPPI.mat']);
    else
        SPM=spm_spm(SPM);
    end
else
    SPM=SPMPPI.SPM;
end

