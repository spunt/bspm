function files = ipctb_ica_rename_4d_file(files, fileNumber)
%For nifti and analyze files rename the files by adding a number at the
%end

if ~isempty(files)

    newFiles = repmat(struct('name', []), size(files, 1), 1);
    for nF = 1:size(files, 1)

        currentFile = deblank(files(nF, :));
        %[pathstr, fileN, extn] = fileparts(currentFile);
        extnPos = ipctb_ica_findstr(currentFile, '.');
        if ~isempty(extnPos)
            extn = currentFile(extnPos(end):end);
            fileN = currentFile(1:extnPos(end)-1);
        else
            fileN = currentFile;
            extn = '';
        end

        if strcmpi(extn, '.nii') || strcmpi(extn, '.img')
            %             [hdr] = ipctb_ica_read_hdr(currentFile);
            %             numFiles = hdr.dime.dim(5);
            try
                ni = ipctb_nifti(currentFile);
                numFiles = ni.dat.dim(4);
            catch
                numFiles = 1;
            end
            if (~exist('fileNumber', 'var')) || (isempty(fileNumber))
                %fileNumber = (1:numFiles);
                tempFileNum = (1:numFiles);
            else
                tempFileNum = fileNumber;
            end

            tempFiles = repmat(struct('name', []), length(tempFileNum), 1);
            for nn = tempFileNum
                if nn <= numFiles
                    tempFiles(nn).name = [currentFile, ',', num2str(nn)];
                else
                    newFiles(nF).name = '';
                end
            end
            newFiles(nF).name = str2mat(tempFiles.name);
            clear tempFiles;
        else
            newFiles(nF).name = currentFile;
        end

    end

    files = cellstr(str2mat(newFiles.name));
    ind = ipctb_ica_good_cells(files);
    files = files(ind);

    if isempty(files)
        files = [];
    else
        files = str2mat(files);
    end

end