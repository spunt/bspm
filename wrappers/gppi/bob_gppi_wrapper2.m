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
gPPIpath = fullfile(getenv('HOME'), 'Github', 'bspm', 'thirdparty', 'gPPI_Toolbox'); % path for gPPI toolbox
studypath =  '/Users/bobspunt/Documents/fmri/lois'; % path for study directory
maskpath = [studypath filesep '_rois_'];
level1path = 'analysis/lois_rwls'; % path of level analysis (relative to subject folder)


% -----------Patterns for finding files and folders -----------
subjectpattern = 'RA*'; % pattern for finding subject folders (use wildcards)
subTAG = 'all'; %  do which subjects? 'all' to do all, position indices to do subset, e.g., [1 3 7]
maskpattern = '*img'; % pattern for finding image files in ROI directory


% -----------General PPI analysis parameters-----------
spmversion = 8; % version of spm (5 or 8)
estimateTAG = 1; % estimate the level 1 models? (1 = yes, 0 = no)
ppi_folder_affix = 'cond'; % string to affix to name of your PPI analysis folder


% -----------Defining the conditions  (psychological effects)-----------
method = 'cond'; % 'trad' for traditional method, 'cond' for generalized (condition-specific) method
conditions ={'High_NS' 'High_Face' 'High_Hand' 'Low_NS' 'Low_Face' 'Low_Hand'}; % conditions to compute PPIs for (must match names used in level 1 analysis)
% {'High_NS' 'High_Face' 'High_Hand' 'Low_NS' 'Low_Face' 'Low_Hand'};
% {'How_Video' 'How_Text' 'Why_Video' 'Why_Text'};
weights = [1 -1]; % weights to apply to conditions in traditional PPI (ignore if doing condition-spec)


% -----------Definining the seed timecourse (physiological effect)-----------
fcontrast = 7; % F contrast to adjust for (corresponds to number of the ess image in level 1 folder)
extract='eig'; % extract timecourse as first eigenvariate ('eig') or mean eigenvariate ('mean')
masks={'mask.img'}; % name of images to use to constrain definition of the seed region
threshold = .5; % values to threshold the images above at (number of values must match number of masks)
VOImin = 10;  % minimum number of voxels to accept a seed region as valid

%-----------Contrasts-----------
convec = [0 1 0 0 0 0 0 0 0 0]; % vector of weights for computing PPI contrast (trad PPI only)
% IMPORTANT! If you are doing generalized PPI, specify your contrasts below

%---------------------------------------------------------------------%
addpath(gPPIpath)

% Find subjects
%---------------------------------------------------------------------%
subdir = files([studypath filesep subjectpattern filesep level1path]);
subnam = regexprep(subdir,level1path,'');
subnam = regexprep(subnam,studypath,'');
subnam = regexprep(subnam,'/','');
nsubs = length(subdir);
if strcmp(subTAG,'all')
    dosubs = 1:nsubs;
else
    dosubs = subTAG;
end

%---------------------------------------------------------------------%
% Begin looping through regions and subjects
%---------------------------------------------------------------------%
for i=dosubs
    cd(gPPIpath)
    load([gPPIpath filesep 'parameters.mat'])
%             P = rmfield(P,'VOI');
    P.subject=subnam{i};
    tmp = files([subdir{i} filesep 'roi_dmpfc*nii']);
    P.VOI.VOI = tmp{1};
    P.VOI.masks = masks;
    P.VOI.thresh = threshold;
    P.VOI.min = VOImin;
    P.VOImin = VOImin;
    P.Region = ['roi_dmpfc' '_' ppi_folder_affix];
    P.SPMver = spmversion;
    P.directory= subdir{i};
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
    for c = 1:length(conditions)
        P.Contrasts(c) = P.Contrasts(1);
        P.Contrasts(c).name = conditions{c};
        P.Contrasts(c).left = conditions(c);
        P.Contrasts(c).right = {'none'};
    end
    %--------------------------------------------------------------%
    end
    % now run! 
    PPPI(P)
end


 
 
 
 
