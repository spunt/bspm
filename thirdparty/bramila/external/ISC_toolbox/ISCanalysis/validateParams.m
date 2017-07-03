function handles = validateParams(handles,vField)

Pub = handles.Pub;

switch vField

    case 'subj'
        filesExist = true;
        % check if source fMRI files exist:
        handles.validFlag = true;
        for k1 = 1:size(Pub.subjectSource,1)
            for k2 = 1:size(Pub.subjectSource,2)
                if ~isempty(Pub.subjectSource{k1,k2})
                    if exist(Pub.subjectSource{k1,k2},'file') ~= 2
                        disp(['  File "' Pub.subjectSource{k1,k2} '" not found.'])
                        filesExist = false;
                        handles.validFlag = false;
                    end                    
                end
            end
        end
        if handles.validFlag
            iter = 1;
            notEm = 0;
            for k = 1:size(Pub.subjectSource,1)
                for m = 1:size(Pub.subjectSource,2)
                    if isempty(Pub.subjectSource{k,m}) == 0
                        notEm = 1;
                        fileExt = Pub.subjectSource{k,m}(end-2:end);
                        if iter > 1
                            if strcmp(fileExt,fileExtPrev) == 0
                                disp('File extension mismatch.')
                                handles.validFlag = false;
                                return
                            end
                        end
                        fileExtPrev = fileExt;
                        iter = iter + 1;
                    end
                end
            end
            if notEm == 1
               if strcmp(fileExt,'nii')
                   Pub.fileFormatSubj = 'nii';
               elseif strcmp(fileExt,'mat')
                   Pub.fileFormatSubj = 'mat';
               else
                   disp('File extension mismatch.')
                   handles.validFlag = false;
               end
            else
                disp('No subjects specified.')
                handles.validFlag = false;
            end
        end
        iter = [];
        for mm = 1:size(Pub.subjectSource,1)
            iter(mm) = 1;
            for kk = 1:size(Pub.subjectSource,2)
                if isempty(Pub.subjectSource{mm,kk}) == 0
                    SS{mm,iter(mm)} = Pub.subjectSource{mm,kk};
                    iter(mm) = iter(mm) + 1;
                end
            end
        end
        if isempty(iter) == 0
            if length(unique(iter)) == 1
                Pub.subjectSource = SS;
                %disp('All subject source files found.')
            else
                disp('Each session must have same number of subjects.')
                handles.validFlag = false;
            end
        else
            disp('No subjects specified.')
            handles.validFlag = false;
        end
        
        % Check sizes of the fMRI data sets. They must all be equal.
        if filesExist
            fileFormat = Pub.fileFormatSubj;
            for k = 1:size(Pub.subjectSource,1)
                dSiz = NaN*ones(size(Pub.subjectSource,2),4);
                for m = 1:size(Pub.subjectSource,2)
                    fileName = Pub.subjectSource{k,m};
                    [siz,flag] = getDataSize(fileName,fileFormat);
                    if ~flag
                        disp('All fMRI data sets must be 4-dimensional.')
                        handles.validFlag = false;
                    end
                    dSiz(m,1:4) = siz;
                end
                if sum(sum(diff(dSiz))) ~= 0
                    disp('All fMRI data sets must have exactly the same number of voxels and time-points in each session.')
                    handles.validFlag = false;
                end
            end
        end
        
        
        if handles.validFlag
            disp('Ok!')
        end
        
        handles.Pub = Pub;
        
    case 'template'
        handles.validFlag = true;
        atype = [{'cort'};{'sub'};{'cort'};{'sub'};{'cort'};{'sub'}];
        atlasTh = [0 0 25 25 50 50];
        
        for k = 1:length(atype)
             brainAtlases = [Pub.atlasPath 'HarvardOxford-'...
             atype{k} '-maxprob-thr' num2str(atlasTh(k)) '-' num2str(handles.Priv.voxelSize) 'mm.nii'];
            if exist(Pub.atlasPath,'dir') ~= 7
                disp([Pub.atlasPath ' is not a directory.'])
                handles.validFlag = false;                
            end
            if exist(brainAtlases,'file') ~= 2
                disp(['  Template "' brainAtlases '" not found.'])
                handles.validFlag = false;    
            end
        end
        brainMask = [Pub.maskPath 'MNI152_T1_' ...
        num2str(handles.Priv.voxelSize) 'mm_brain_mask.nii'];
        if exist(Pub.maskPath,'dir') ~= 7
            disp([Pub.maskPath ' is not a directory.'])
            handles.validFlag = false;                    
        end
        if exist(brainMask) ~= 2
            disp(['  Template "' brainMask '" not found.'])
            handles.validFlag = false;
        end
            
        if handles.validFlag == 1
            if handles.validFlag
                disp('Ok!')
            end
        end

    case 'freqBands'
        N = min(handles.Priv.dataSize(:,4));
        lw = 4;
        maxLev = fix(log(N/(lw-1))/log(2)) - 1;
        %maxLev = wmaxlev(handles.Priv.dataSize,'db2') - 1;
        if handles.Pub.nrFreqBands > maxLev
            disp(['Maximum number of frequency bands is ' num2str(maxLev) '.'])
            handles.validFlag = false;
            return
        end
        if handles.validFlag
            disp('Ok!')
        end
    case 'TR'
        if isnan(handles.Pub.samplingFrequency)
            disp('TR of the fMRI time-series must be provided (in seconds)!')
            handles.validFlag = false;
            return 
        end
        
end
