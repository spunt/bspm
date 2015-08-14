function subFolders = ipctb_ica_listSubFolders(sourceDir)
% list sub folders

% list sub folders
subFolders = ipctb_ica_listDir(sourceDir);
tempVec = [];

index1 = strmatch('.', subFolders, 'exact');
index2 = strmatch('..', subFolders, 'exact');

index = [index1 index2];

CheckOne = ones(size(subFolders, 1), 1);

CheckOne(index) = 0;

tempVec = find(CheckOne ~= 0);

if ~isempty(tempVec)
    subFolders = subFolders(tempVec, :);
    for ii = 1:size(subFolders, 1)
        store(ii).name = deblank(fullfile(sourceDir, subFolders(ii, :)));
    end

    clear  subFolders;
    subFolders = str2mat(store.name);
    clear store;
else
    subFolders = [];
end



