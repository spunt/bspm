function out = spm_rwls_run_fmri_est(job)
% Set up the design matrix and run a design.
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_rwls_run_fmri_est.m 3327 2010-05-18 08:27:32Z joern $


global defaults
if isempty(defaults)
    spm_defaults;
end;
if ~isfield(defaults,'modality')
    defaults.modality = 'FMRI';
end;

%-Load SPM.mat file
%-----------------------------------------------------------------------
SPM = [];
load(job.spmmat{:});
% pre-set output
out.spmmat = job.spmmat;

original_dir = pwd;

%-Move to the directory where the SPM.mat file is
%-----------------------------------------------------------------------
cd(fileparts(job.spmmat{:}));



%=======================================================================
% R E M L   E S T I M A T I O N
%=======================================================================
if isfield(job.method,'Classical'),
    
    SPM = spm_rwls_spm(SPM);
    
    %-Automatically set up contrasts for factorial designs
    %-------------------------------------------------------------------
    if isfield(SPM,'factor')
        if SPM.factor(1).levels > 1
        % don't both if you've only got 1 level and 1 factor
            cons = spm_design_contrasts(SPM);
        
            %-Create F-contrasts
            %-----------------------------------------------------------
            for i=1:length(cons)
                con  = cons(i).c;
                name = cons(i).name;
                STAT = 'F';
                [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
                if all(I)
                    DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                else
                    DxCon = [];
                end
                if isempty(SPM.xCon),
                    SPM.xCon = DxCon;
                else
                    SPM.xCon(end+1) = DxCon;
                end
               SPM = spm_contrasts(SPM,length(SPM.xCon));
            end
        
            %-Create t-contrasts
            %-----------------------------------------------------------
            for i=1:length(cons)
                % Create a t-contrast for each row of each F-contrast
                % The resulting contrast image can be used in a 2nd-level analysis
                Fcon  = cons(i).c;
                nrows = size(Fcon,1);
                STAT  = 'T';
                for r=1:nrows,
                    con = Fcon(r,:); 
                    str = cons(i).name;
                    if ~isempty(strmatch('Interaction',str))
                        name = ['Positive ',str,'_',int2str(r)];
                    else
                        sp1  = min(find(isspace(str))); 
                        name = ['Positive',str(sp1:end),'_',int2str(r)];
                    end
                
                    [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
                    if all(I)
                        DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                    else
                        DxCon = [];
                    end
                    if isempty(SPM.xCon),
                        SPM.xCon = DxCon;
                    else
                        SPM.xCon(end+1) = DxCon;
                    end
                   SPM = spm_contrasts(SPM,length(SPM.xCon));
                end
            end
        end % if SPM.factor(1).levels > 1
    end % if isfield(SPM,'factor')
    
    %out.spmvar = SPM;
    out.beta = cellfun(@(fn)fullfile(SPM.swd,fn), cellstr(char(SPM.Vbeta(:).fname)),'UniformOutput',false);
    out.mask = {fullfile(SPM.swd,SPM.VM.fname)};
    out.resms = {fullfile(SPM.swd,SPM.VResMS.fname)};
    cd(original_dir); % Change back
    fprintf('Done\n');
    return
end


