function removeMemmapData(Params)
% The function removes the wanted memmapped binary data files (as defined in the GUI. 
% If the files are deleted the program generates them again when the analysis is rerun. 
%Note that regeneration of the filtered data can take a long time.

%Priv = Params.PrivateParams;
Pub = Params.PublicParams;

if(Pub.removeMemmaps)
    disp('Removing memmapped source data:')
    %get the current path ending (/ or \)
    slash = Pub.dataDestination(end);
    %define path and find the bin files
    rmpath = [Pub.dataDestination, 'fMRIpreprocessed', slash];
    files = dir([rmpath, '*.bin']);
    %deleting files one by one, just in case
    for k = 1:length(files)
        disp(['Removing file: ' files(k).name])
        delete([rmpath, files(k).name])
        fprintf('\b',1); disp('... Done!')
    end
    %should these also get cleared from memmap struct? 
else
    disp('Memmap removal not selected, skipping removal phase')
end

if(Pub.removeFiltermaps)
    disp('Removing Filtered data:')
    %get the current path ending (/ or \)
    slash = Pub.dataDestination(end);
    %define path and find the bin files
    rmpath = [Pub.dataDestination, 'fMRIfiltered', slash];
    files = dir([rmpath, '*.bin']);
    %deleting files one by one, just in case
    for k = 1:length(files)
        disp(['Removing file: ' files(k).name])
        delete([rmpath, files(k).name])
        fprintf('\b',1); disp('... Done!')
    end
    %should these also get cleared from memmap struct? 
else
    disp('Filtermap removal not selected, skipping removal phase')
end