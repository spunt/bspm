function [memMaps Params] = changeDirectory(parentDir,atlasPath)

% This function changes file pointers in all memory-map objects and changes
% all necessary parameter -fields in parameter-struct. Function is a part of
% ISC toolbox and must be run if analysis data has been moved from its original
% location where memory mapping was performed.
%
% inputs:
%
% parentDir - full directory path name where analysis data is located.
%             e.g. if new data location is:
%                 /home/data/
%                 /home/data/results
%                 /home/data/fMRIpreprocessed
%                 /home/data/fMRIfiltered
%                 .......................
%
%             You must set: parentDir = /home/data/
%
% atlasPath - full directory path where (FSL nifti-format) atlas data is located
%
%
% outputs:
% memMaps - struct containing updated memory maps (data can be accesses through them)
% Params - updated analysis parameter -struct
%
% Note!! Outputs are optional and are needed just if user wants to access data and
% investigate Parameter changes via Matlab's workspace. Updated variables are automatically
% saved in the parentDir so there is no need to save the variables after running this function.

% Jukka-Pekka Kauppi
% 05.09.2010



if nargin > 2
    error('Number of inputs must be 2!!')
    return
else
    if nargin == 1
        atlasPath = '';
    end
    [flag,directories,Params,memMaps] = checkInputs(parentDir,atlasPath);
    if ~flag
        return
    end
end

maskPath = atlasPath;

Priv = Params.PrivateParams;
Pub = Params.PublicParams;
disp(' ')
disp('Change mask path.....')

try
    Priv.brainMask = [maskPath 'MNI152_T1_' ...
        num2str(Priv.voxelSize) 'mm_brain_mask.nii'];
    atype = [{'cort'};{'sub'};{'cort'};{'sub'};{'cort'};{'sub'}];
    Pub.maskPath = maskPath;
    disp('OK!')
catch
    disp('Path not specified or it''s not found, update ignored!')
end
disp(' ')
disp('Change brain atlas path.....')
try
    for k = 1:length(atype)
        Priv.brainAtlases{k} = [atlasPath 'HarvardOxford-'...
            atype{k} '-maxprob-thr' num2str(Priv.atlasTh(k)) ...
            '-' num2str(Priv.voxelSize) 'mm.nii'];
    end
    Pub.atlasPath = atlasPath;
    disp('OK!')
catch
    disp('Path not specified or it''s not found, update ignored!')
end

Pub.dataDestination = parentDir;
Priv.PFDestination = [parentDir 'PF' parentDir(end)];
Priv.statsDestination = [parentDir 'stats' parentDir(end)];
Priv.subjectDestination = [parentDir 'fMRIpreprocessed' parentDir(end)];
Priv.subjectFiltDestination = [parentDir 'fMRIfiltered' parentDir(end)];
Priv.resultsDestination = [parentDir 'results' parentDir(end)];

R(1) = Pub.ssiOn;
R(2) = Pub.nmiOn;
R(3) = Pub.corOn;
R(4) = Pub.kenOn;
R = nonzeros(R.*(1:4));

endS = [{'.bin'},{'_win.bin'}];
disp(' ')
disp('Update memory map pointers.....')

disp(' ')
disp('Synchronization maps:')
try
    % update synchronization data parentDir names:
    fn = {'whole','win'};
    for k = 1:length(fn)
        for m = 0:Priv.maxScale + 1
            for n = 1:Priv.nrSessions
                for p = 1:length(R)
                    % get memMap:
                    H = memMaps.(Priv.resultMapName).(fn{k}).([Priv.prefixFreqBand ...
                        num2str(m)]).([Priv.prefixSession num2str(n)]).(Priv.simM{R(p)});
                    
                    % set parentDirName:
                    H.fileName = [parentDir directories{1} Priv.prefixResults '_' ...
                        Priv.simM{R(p)} '_' Priv.prefixFreqBand ...
                        num2str(m) '_' Priv.prefixSession num2str(n) '_' ...
                        Priv.transformType endS{k}];
                    
                    memMaps.(Priv.resultMapName).(fn{k}).([Priv.prefixFreqBand ...
                        num2str(m)]).([Priv.prefixSession ...
                        num2str(n)]).(Priv.simM{R(p)}) = H;
                end
            end
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Maps not found, update ignored!')
end

disp(' ')
disp('Preprocessed data:')
try
    % update original data matrix parentDir names:
    for n = 1:Priv.nrSubjects
        for m = 1:Priv.nrSessions
            % get memMap:
            H2 = memMaps.(Priv.origMapName).([Priv.prefixSession ...
                num2str(m)]).([Priv.prefixSubject num2str(n)]);
            
            % set parentDirName:
            H2.fileName = [parentDir directories{3} Priv.prefixSubject num2str(n) ...
                Priv.prefixSession num2str(m) endS{1}];
            
            memMaps.(Priv.origMapName).([Priv.prefixSession ...
                num2str(m)]).([Priv.prefixSubject num2str(n)]) = H2;
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Data not found, update ignored!')
end


disp(' ')
disp('Filtered data:')
try
    % update filtered data matrix parentDir names:
    for n = 1:Priv.nrSubjects
        for m = 1:Priv.nrSessions
            for p = 1:Priv.maxScale + 1
                % get memMap:
                H3 = memMaps.(Priv.filtMapName).([Priv.prefixSession ...
                    num2str(m)]).([Priv.prefixSubjectFilt ...
                    num2str(n)]).([Priv.prefixFreqBand num2str(p)]);
                
                % set parentDir name:
                H3.fileName = [parentDir directories{2} ...
                    Priv.prefixSubjectFilt num2str(n) ...
                    '_' Priv.prefixFreqBand num2str(p) '_' ...
                    Priv.prefixSession num2str(m) '_' ...
                    Priv.transformType '.bin'];
                
                memMaps.(Priv.filtMapName).([...
                    Priv.prefixSession ...
                    num2str(m)]).([Priv.prefixSubjectFilt ...
                    num2str(n)]).([Priv.prefixFreqBand ...
                    num2str(p)]) = H3;
            end
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Data not found, update ignored!')
end

disp(' ')
disp('Phase maps:')
try
    % update phase map parentDir names:
    for m = 1:Priv.nrSessions
        for p = 0:Priv.maxScale + 1
            % get memMap:
            H3 = memMaps.(Priv.phaseMapName).([Priv.prefixSession ...
                num2str(m)]).([Priv.prefixFreqBand num2str(p)]);
            
            % set parentDir name:
            H3.fileName = [parentDir directories{6} ...
                'phase_' Priv.prefixFreqBand num2str(p) '_' ...
                Priv.prefixSession num2str(m) '_' ...
                Priv.transformType '.bin'];
            
            memMaps.(Priv.phaseMapName).([Priv.prefixSession ...
                num2str(m)]).([Priv.prefixFreqBand num2str(p)]) = H3;
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Data not found, update ignored!')
end

disp(' ')
disp('Intersubject synchronization curves:')
try
    % update synchronization curve parentDir names:
    for n = 1:Priv.nrSessions
        for m = 0:Priv.maxScale + 1
            memMaps.(Priv.synchMapName).(...
                [Priv.prefixSession num2str(n)]).([Priv.prefixFreqBand ...
                num2str(m)]).fileName = [parentDir directories{1} Priv.prefixSyncResults ...
                Priv.prefixSession num2str(n) ...
                Priv.prefixFreqBand num2str(m) '.bin'];
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Curves not found, update ignored!')
end

disp(' ')
disp('Phase synchronization curves:')
try
    % update synchronization curve parentDir names:
    for n = 1:Priv.nrSessions
        for m = 0:Priv.maxScale + 1
            memMaps.(Priv.phaseSynchMapName).(...
                [Priv.prefixSession num2str(n)]).([Priv.prefixFreqBand ...
                num2str(m)]).fileName = [parentDir directories{6} Priv.prefixPhaseSyncResults ...
                Priv.prefixSession num2str(n) ...
                Priv.prefixFreqBand num2str(m) '.bin'];
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Curves not found, update ignored!')
end



disp(' ')
disp('ZPF maps:')
try
    
    % update ZPF parentDir names:
    fc = ((Priv.maxScale+2)^2-(Priv.maxScale+2))/2;
    %fn = {'whole','win'};
    fn = {'whole'};
    for k = 1:length(fn)
        for n = 1:Priv.nrSessions
            for p = 1:fc
                % get memMap:
                H = memMaps.(Priv.PFMapName).(fn{k}).([...
                    Priv.prefixSession num2str(n)]).(...
                    Priv.simM{3}).([Priv.prefixFreqComp num2str(p)]);
                % set parentDirName:
                H.fileName = [parentDir directories{5} Priv.prefixPF '_' ...
                    Priv.simM{3} '_' Priv.prefixSession num2str(n) '_' ...
                    Priv.transformType Priv.prefixFreqComp num2str(p) endS{k}];
                memMaps.(Priv.PFMapName).(fn{k}).([...
                    Priv.prefixSession num2str(n)]).(...
                    Priv.simM{3}).([Priv.prefixFreqComp num2str(p)]) = H;
            end
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Maps not found, update ignored!')
end

disp(' ')
disp('Stats maps:')
try
    % update stats data parentDir names:
    fn = {'whole','win'};
    for k = 1:length(fn)
        for m = 0:Priv.maxScale + 1
            for n = 1:Priv.nrSessions
                % get memMap:
                H = memMaps.(Priv.statMapName).(fn{k}).([Priv.prefixFreqBand ...
                    num2str(m)]).([Priv.prefixSession num2str(n)]).(Priv.simM{3});
                % set parentDirName:
                H.fileName = [parentDir directories{4} Priv.prefixTMap '_' ...
                    Priv.simM{3} '_' Priv.prefixFreqBand ...
                    num2str(m) '_' Priv.prefixSession num2str(n) '_' ...
                    Priv.transformType endS{k}];
                memMaps.(Priv.statMapName).(fn{k}).([Priv.prefixFreqBand ...
                    num2str(m)]).([Priv.prefixSession num2str(n)]).(Priv.simM{3}) = H;
            end
        end
    end
    disp('OK!')
catch
    disp(lasterr)
    disp('Maps not found, update ignored!')
end
disp(' ')
disp('Correlation matrices:')
try
    % update correlation matrix data parentDir names:
    fn = {'whole','win'};
    for k = 1:length(fn)
        for m = 0:Priv.maxScale + 1
            for n = 1:Priv.nrSessions
                % get memMap:
                H = memMaps.(Priv.cormatMapName).(fn{k}).([Priv.prefixFreqBand ...
                    num2str(m)]).([Priv.prefixSession num2str(n)]).(Priv.simM{3});
                % set parentDirName:
                H.fileName = [parentDir directories{4} Priv.prefixCorMat '_' ...
                    Priv.simM{3} '_' Priv.prefixFreqBand ...
                    num2str(m) '_' Priv.prefixSession num2str(n) '_' ...
                    Priv.transformType endS{k}];
                memMaps.(Priv.cormatMapName).(fn{k}).([Priv.prefixFreqBand ...
                    num2str(m)]).([Priv.prefixSession num2str(n)]).(Priv.simM{3}) = H;
            end
        end
    end
    
    
    
    disp('OK!')
catch
    disp(lasterr)
    disp('Matrices not found, update ignored!')
end



disp(' ')
Params.PublicParams = Pub;
Params.PrivateParams = Priv;
disp(['Saving updated variables to ' parentDir])

save([Pub.dataDestination Pub.dataDescription '.mat'],'Params')
save([Pub.dataDestination 'memMaps.mat'],'memMaps')

disp(' ')
disp('done!')





function [flag,directories,Params,memMaps] = checkInputs(parentDir,atlasPath)
disp('Checking inputs.....')
flag = true;
if ( ~ischar(parentDir) || ~ischar(atlasPath) )
    error('Inputs must be strings!!')
    flag = false;
end
if ~strcmp(parentDir(end),'/') && ~strcmp(parentDir(end),'\')
    error('Parent directory name must end with slash!!')
    flag = false;
    return
end
if ~strcmp(atlasPath(end),'/') && ~strcmp(atlasPath(end),'\')
    error('Atlas directory name must end with slash!!')
    flag = false;
    return
end

if exist(parentDir,'dir') ~= 7
    error('Parent directory does not exist!!')
    flag = false;
    return
end
directories = {['results' parentDir(end)],['fMRIfiltered' parentDir(end)],...
    ['fMRIpreprocessed' parentDir(end)],['stats' parentDir(end)],...
    ['PF' parentDir(end)],['phase' parentDir(end)]};
if exist([parentDir directories{1}],'dir') ~= 7
    mkdir([parentDir directories{1}])
end
if exist([parentDir directories{2}],'dir') ~= 7
    mkdir([parentDir directories{2}])
end
if exist([parentDir directories{3}],'dir') ~= 7
    mkdir([parentDir directories{3}])
end
if exist([parentDir directories{4}],'dir') ~= 7
    mkdir([parentDir directories{4}])
end
if exist([parentDir directories{5}],'dir') ~= 7
    mkdir([parentDir directories{5}])
end
if exist([atlasPath],'dir') ~= 7
    disp('Invalid atlas path!!')
end

disp('Searching Params-struct from parent directory.....')
s = what(parentDir);
fl = false;
for k = 1:length(s.mat)
    q = whos('-file',[parentDir s.mat{k}]);
    for h = 1:length(q)
        if strcmp(q(h).name,'Params');
            load([parentDir s.mat{k}],q(h).name)
            if isfield(Params,'PublicParams')
                if isfield(Params.PublicParams,'dataDescription')
                    if strcmp([Params.PublicParams.dataDescription '.mat'],s.mat{k});
                        disp(['  Found ' s.mat{k}])
                        fl = true;
                        break
                    end
                end
            end
        end
    end
    if fl
        break
    end
    
end
if ~fl
    error('mat-file containing parameter-struct not found in parent directory!!')
    flag = false;
end
disp('Searching memMap.m from parent directory.....')
fl = false;
for k = 1:length(s.mat)
    if strcmp(s.mat{k},'memMaps.mat')
        load([parentDir s.mat{k}])
        fl = true;
    end
end
if ~fl
    error('memMaps.mat not found in parent directory!!')
    flag = false;
end
disp('  Found memMaps.mat')
