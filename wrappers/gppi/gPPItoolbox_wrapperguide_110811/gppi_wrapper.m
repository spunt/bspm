%================================================%
% GPPI_WRAPPER
%
%   This script allows you to run Donald McLaren's gPPI toolbox to build 
%   and estimate psychophysiological interaction (PPI) models for multiple
%   subjects and for multiple seed volumes-of-interest (VOI). Note that
%   before using this, you must first create image files corresponding to
%   your seed regions. In addition, note that there are TWO sections which
%   require your input. In the first section, you will define the majority
%   of the variables necessary for performing the analysis. In the second
%   section (near the bottom of the script), you will need to specify the
%   contrasts that you would like to estimate. 
%   
%   You can download the toolbox and find further information at:
%   http://www.brainmap.wisc.edu/PPI
%
%   November 7, 2011 -- Created by Bob Spunt
%================================================%
clear all; home;

% USER INPUT: Paths, file-finding patterns, and analysis parameters
%---------------------------------------------------------------------%

% -----------Paths-----------
gPPIpath = '/Applications/fmri/data/matchrest/_scripts_/gPPItoolbox'; % path for gPPI toolbox
studypath = '/Applications/fmri/data/matchrest'; % path for study directory
maskpath = '/Applications/fmri/data/matchrest/_rois_'; % path for folder containing ROIs
level1path = 'analysis/matchrest'; % path of level analysis (relative to subject folder)


% -----------Patterns for finding files and folders -----------
subjectpattern = 'MBC*'; % pattern for finding subject folders (use wildcards)
subTAG = 'all'; %  do which subjects? 'all' to do all, position indices to do subset, e.g., [1 3 7]
maskpattern = '*img'; % pattern for finding image files in ROI directory


% -----------General PPI analysis parameters-----------
spmversion = 8; % version of spm (5 or 8)
estimateTAG = 1; % estimate the level 1 models? (1 = yes, 0 = no)
ppi_folder_affix = 'cond_110811'; % string to affix to name of your PPI analysis folder


% -----------Defining the conditions  (psychological effects)-----------
method = 'cond'; % 'trad' for traditional method, 'cond' for generalized (condition-specific) method
conditions = {'Rest' 'Match'}; % conditions to compute PPIs for (must match names used in level 1 analysis)
weights = [1 -1]; % weights to apply to conditions in traditional PPI (ignore if doing condition-spec)


% -----------Definining the seed timecourse (physiological effect)-----------
fcontrast = 2; % F contrast to adjust for (corresponds to number of the ess image in level 1 folder)
extract='eig'; % extract timecourse as first eigenvariate ('eig') or mean eigenvariate ('mean')
masks={'mask.img'}; % name of images to use to constrain definition of the seed region
threshold = [.5]; % values to threshold the images above at (number of values must match number of masks)
VOImin = 1;  % minimum number of voxels to accept a seed region as valid


%-----------Contrasts-----------
convec = [0 1 0 0 0 0 0 0 0 0]; % vector of weights for computing PPI contrast (trad PPI only)
% IMPORTANT! If you are doing generalized PPI, specify your contrasts below

%---------------------------------------------------------------------%
addpath(gPPIpath)

% Find subjects
%---------------------------------------------------------------------%
cd(studypath);
fprintf('\nSUBJECT LIST:\n');
d=dir(subjectpattern);
for i=1:length(d)
    subnam{i}=d(i).name;
    subdir{i} = [studypath '/' subnam{i}];
    fprintf('\tAdding %s to subject list\n',subnam{i})
end
nsubs = length(subnam);
if strcmp(subTAG,'all')
    dosubs = 1:nsubs;
else
    dosubs = subTAG;
end

% Find mask files
%---------------------------------------------------------------------%
cd(maskpath)
fprintf('\n\nMASK LIST:\n');
d=dir(maskpattern);
for i=1:length(d)
    voiNAMES{i}=d(i).name;
    voi_anal_names{i} = voiNAMES{i}(1:end-4);
    fprintf('\tAdding %s to mask list\n',voiNAMES{i})
end
nmasks = length(voiNAMES);
domasks = 1:nmasks;
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
% Begin looping through regions and subjects
%---------------------------------------------------------------------%
for v=domasks
    for i=dosubs
        try
            cd(gPPIpath)
            load('parameters.mat')
            P.subject=subnam{i};
            P.VOI.VOI = [maskpath filesep voiNAMES{v}];
            P.VOI.masks = masks;
            P.VOI.thresh = threshold;
            P.VOI.min = VOImin;
            P.VOImin = VOImin;
            P.Region = [voi_anal_names{v} '_' ppi_folder_affix];
            P.SPMver = spmversion;
            P.directory=[subdir{i} filesep level1path];
            P.method = method;
            P.Estimate = estimateTAG;
            P.contrast = fcontrast;
            P.extract = extract;     
            if strcmp(method,'trad')
                P.Tasks = conditions;
                P.Weights = weights;
                P.Contrasts(1).name = 'PPI';
                P.Contrasts(1).left = convec;
                P.Contrasts(1).right = [];
            else
                P.Tasks = ['0' conditions];
            end
            
            if strcmp(method,'cond')
                
            % USER INPUT: Specify your contrasts of interest (Cond PPI only)
            %--------------------------------------------------------------%
            
                P.Contrasts(1).name = 'Rest-Match';
                P.Contrasts(1).left = {'Rest'};
                P.Contrasts(1).right = {'Match'};

                P.Contrasts(2) = P.Contrasts(1);
                P.Contrasts(2).name = 'Rest-Base';
                P.Contrasts(2).left = {'Rest'};
                P.Contrasts(2).right = {''};    

                P.Contrasts(3) = P.Contrasts(1);
                P.Contrasts(3).name = 'Match-Base';
                P.Contrasts(3).left = {'Match'};
                P.Contrasts(3).right = {''};
                
%                 P.Contrasts(4) = P.Contrasts(1);
%                 P.Contrasts(4).name = 'HowText-Base';
%                 P.Contrasts(4).left = {'HowVerbal'};
%                 P.Contrasts(4).right = {''};

            %--------------------------------------------------------------%
            end
            % now run! 
            PPPI(P)
            
        catch
            disp(['Failed: ' subnam{i}])
        end
        clear P
    end
end

