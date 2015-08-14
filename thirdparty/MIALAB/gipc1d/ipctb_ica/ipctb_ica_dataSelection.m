function [files, designMatrix, numOfSub, numOfSess, dataSelMethod, diffTimePoints, spmMatFlag]  = ipctb_ica_dataSelection(...
    inputFile, outputDir, output_prefix, read_complex_file_naming, read_complex_images)
% Data selection step in GIFT toolbox
%
% Inputs:
% 1. inputFile - inputFile containing the directory and file pattern for
% subjects
% 2. outputDir - Output directory for the analysis.
% 3. output_prefix - File output prefix
% 4. read_complex_file_naming - File naming for reading complex images
% 5. read_complex_images - Complex data type (real&imaginary) or
% (magnitude&phase)
%
% Output:
% 1. files - data structure containing the subject images. Data will be stored
% in sessions.
% 2. designMatrix - design matrix structure
% 3. numOfSub - Number of subjects selected
% 4. numOfSess - Number of sessions selected
% 5. dataSelMethod - Data selection method type
% 6. diffTimePoints - Time points vector.
% 7. spmMatFlag - flag to distinguish how to select the design matrix.

ipctb_ica_defaults;
global FUNCTIONAL_DATA_FILTER;

if ~exist('inputFile', 'var')
    inputFile = [];
end

% Get modality type
[modalityType, dataTitle] = ipctb_ica_get_modality;


% Initialise output vars
files.name = [];
designMatrix.name = [];
dataSelMethod = 2;
diffTimePoints = [];
spmMatFlag = 'no';
numOfSub = 0;
numOfSess = 0;

if ~exist('outputDir', 'var')
    outputDir = pwd;
end

if ~exist('output_prefix', 'var')
    output_prefix = '';
end

if ~exist('inputFile', 'var')
    inputFile = [];
end

if ~exist('dataType', 'var')
    dataType = 'real';
end


fileNumbers = [];

% get the subject matrix file naming
subjectFile = [output_prefix, 'Subject.mat'];
subjectFile = fullfile(outputDir, subjectFile);

[modalityType, dataTitle] = ipctb_ica_get_modality;

if isempty(inputFile)
    % Use GUI to select the data selection method

    % open a popup window asking the user to select the data
    questionString = 'Is your data stored in one folder?';
    choiceString = str2mat('Yes', 'No');
    datasel_Handle = ipctb_ica_getGraphics('Selecting data method', 'normal', 'data selection'); % figure handle
    set(datasel_Handle, 'menubar', 'none');


    % plot menu on this figure window
    % Help Menu
    helpMenu = uimenu('parent', datasel_Handle, 'label', 'GIFT-Help');
    htmlHelpMenu = uimenu(helpMenu, 'label', 'Data Selection', 'callback', ...
        'ipctb_ica_openHTMLHelpFile(''screenshots_example_data.htm'');');

    popupAnswer = ipctb_ica_promptUI('popup', questionString, choiceString, 'numeric', datasel_Handle);
    delete(datasel_Handle);

    % Data selection method
    if popupAnswer == 1
        dataSelMethod = 1;
    else
        dataSelMethod = 2;
    end

else

    % Read the data selection method
    keywd = 'dataSelectionMethod';
    try
        inputData = ipctb_ica_read_variables(inputFile, keywd, 'scalar', 'integer');
        % get the data selection method field
        dataSelMethod = getfield(inputData, keywd);
        clear inputData;

        % Validate data selection method
        if (dataSelMethod ~= 1) & (dataSelMethod ~= 2)
            dataSelMethod = 2;
            disp('dataSelMethod variable takes only values 1 and 2. By default selecting value 2...');
        end
        % set the input text for controls
    catch
    end

    % design matrix
    keywd = 'keyword_designMatrix';
    inputData = ipctb_ica_read_variables(inputFile, keywd, 'scalar', 'character');
    spmMatFlag = lower(getfield(inputData, keywd));
    clear inputData;

    % change the spm mat flag
    if strcmpi(spmMatFlag, 'one')
        spmMatFlag = 'same_sub_same_sess';
    end

    % change the spm mat flag
    if strcmpi(spmMatFlag, 'all')
        spmMatFlag = 'diff_sub_diff_sess';
    end

    % get the spm mat flag
    if ~strcmp(spmMatFlag, 'same_sub_same_sess') & ~strcmp(spmMatFlag, 'same_sub_diff_sess') & ...
            ~strcmp(spmMatFlag, 'diff_sub_diff_sess') & ~strcmp(spmMatFlag, 'no')
        error(['Please check the keyword_designMatrix variable in input file: ', inputFile]);
    end

    % For all subjects search for the suffix designMat
    if  strcmpi(spmMatFlag, 'same_sub_diff_sess') | strcmpi(spmMatFlag, 'same_sub_same_sess')
        keywd = 'OnedesignMat';
        inputData = ipctb_ica_read_variables(inputFile, keywd, 'scalar', 'file');
        designMatrix.name = getfield(inputData, keywd);
        clear inputData;
    end

end


% answer to question 1
if dataSelMethod == 1

    if isempty(inputFile)
        data_setDir = ipctb_ica_selectEntry('typeEntity', 'directory', 'title', ...
            'Select root folder for subjects and sessions');
        if ~isempty(data_setDir)
            % dialog Title
            dlg_title = [dataTitle, ' data information.'];

            numParameters = 1;

            %                 inputText(numParameters).promptString = 'Select data type.';
            %                 inputText(numParameters).uiType = 'popup';
            %                 inputText(numParameters).answerString = {'Real', 'Complex'};
            %                 inputText(numParameters).dataType = 'string';
            %                 inputText(numParameters).tag = 'datatype';
            %                 inputText(numParameters).enable = 'on';
            %
            %                 numParameters = numParameters + 1;

            inputText(numParameters).promptString = 'Select file pattern for reading data.';
            inputText(numParameters).uiType = 'edit';
            inputText(numParameters).answerString = FUNCTIONAL_DATA_FILTER;
            inputText(numParameters).dataType = 'string';
            inputText(numParameters).tag = 'filepattern';
            inputText(numParameters).enable = 'on';

            numParameters = numParameters + 1;

            inputText(numParameters).promptString = 'Are session folders inside subject folders?';
            inputText(numParameters).uiType = 'popup';
            inputText(numParameters).answerString = str2mat('Yes', 'No');
            inputText(numParameters).dataType = 'string';
            inputText(numParameters).tag = 'data_folder';
            inputText(numParameters).enable = 'on';

            numParameters = numParameters + 1;

            inputText(numParameters).promptString = 'Enter file numbers to include. Leave empty if you want to select all.';
            inputText(numParameters).uiType = 'edit';
            inputText(numParameters).answerString = '';
            inputText(numParameters).dataType = 'string';
            inputText(numParameters).tag = 'file_numbers';
            inputText(numParameters).enable = 'on';

            %numParameters = numParameters + 1;

            %                 inputText(numParameters).promptString = 'Select an option to read complex images.';
            %                 inputText(numParameters).uiType = 'popup';
            %                 inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
            %                 inputText(numParameters).dataType = 'string';
            %                 inputText(numParameters).tag = 'read_complex_images';
            %                 inputText(numParameters).enable = 'off';
            %
            %                 numParameters = numParameters + 1;
            %
            %
            %                 inputText(numParameters).promptString = 'Select an option to write complex images.';
            %                 inputText(numParameters).uiType = 'popup';
            %                 inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
            %                 inputText(numParameters).dataType = 'string';
            %                 inputText(numParameters).tag = 'write_complex_images';
            %                 inputText(numParameters).enable = 'off';

            numUIControls = length(inputText);

            % Input dialog box (get the necessary numbers)
            answer = ipctb_ica_inputDialog('inputtext', inputText, 'Title', dlg_title);

            getSPMMatrix = 'no';
            if ~isempty(answer)
                % data type
                %dataType = answer{1};
                % file pattern
                filePattern = answer{1};
                if strcmpi(answer{2}, 'yes')
                    % flag folder
                    flagFolder = 'data_in_subject_subfolder';
                else
                    flagFolder = 'data_in_subject_folder';
                end
                % read complex images
                %read_complex_images = answer{4};
                % write complex images
                %write_complex_images = answer{5};
            else
                error('information regarding data type and images should be specified');
            end

        else
            error('data-sets directory is not selected');
        end

        designMatrix.name = [];

        try
            if ~isempty(answer{3})
                fileNumbers = str2num(answer{3});
            end
        catch
            disp('File numbers are not entered correctly');
        end

    else

        % read data-sets directory and filter file pattern
        % keyword for entering subject data
        keywd = 'sourceDir_filePattern_flagLocation';
        inputData = ipctb_ica_read_variables(inputFile, keywd, 'vector', 'cell');
        sourceDir_fileP_flag = getfield(inputData, keywd);
        clear inputData;
        data_setDir = sourceDir_fileP_flag{1};
        filePattern = sourceDir_fileP_flag{2};
        flagFolder = sourceDir_fileP_flag{3};

        if length(sourceDir_fileP_flag) == 4
            fileNumbers = sourceDir_fileP_flag{4};
        end

    end
    % end for getting the information through input file or GUI

    drawnow;
    disp(['Reading data from source directory ', data_setDir]);
    % get the data_sets automatically
    [files, numOfSub, numOfSess, selected_data_sets, data_folders] = ipctb_ica_get_sub_data(data_setDir, filePattern, ...
        flagFolder);

    % There are two ways to select the data get design matrix for each subject
    if strcmpi(spmMatFlag, 'diff_sub_diff_sess')

        keywd = 'spmDesignFilter';
        inputData = ipctb_ica_read_variables(inputFile, keywd, 'scalar', 'string');
        spmDesignFilter = getfield(inputData, keywd);
        clear inputData;

        % loop over subjects
        for ii = 1:numOfSub
            % current subject folder
            currentSub = deblank(data_folders((ii-1)*numOfSess + 1).name);
            % design matrix
            tempdesignMat = ipctb_ica_listFiles_inDir(currentSub, spmDesignFilter);
            if isempty(tempdesignMat)
                try
                    tempdesignMat = ipctb_ica_get_sub_data(currentSub, spmDesignFilter, 'data_in_subject_folder', 0);
                catch
                    error(['Cannot find SPM design matrix for subject ', num2str(ii), ...
                        ' in its folder or sub-folders where folder path is ', ...
                        currentSub]);
                end

                % get the full file path for SPM design matrix
                designMatrix(ii).name = deblank(tempdesignMat(1).name(1, :));
            else
                % get the full file path for SPM design matrix
                designMatrix(ii).name = fullfile(currentSub, deblank(tempdesignMat(1, :)));
            end
            clear tempdesignMat;
        end
        % end loop over subjects
    end
    % end for getting the design matrix

    selectedSubTxtFile = fullfile(outputDir, [output_prefix, 'SelectedDataFolders.txt']);
    % open the file for recording data-sets
    fid = fopen(selectedSubTxtFile, 'w');
    for nLines = 1:size(selected_data_sets, 1)
        tempStr = deblank(selected_data_sets(nLines, :));
        fprintf(fid, '%s\n', tempStr);
    end
    fclose(fid);
    fprintf('\n');
    % closing the file
    disp(['Please see the text file ', selectedSubTxtFile, ' for the selected data folders in order']);

    nF = 0;
    % Loop over subjects
    for nSub = 1:numOfSub
        % Loop over sessions
        for nSess = 1:numOfSess
            % Loop over files
            % for nF = 1:length(files)
            nF = nF + 1;
            tempFiles = ipctb_ica_rename_4d_file(files(nF).name);
            currentFileNum = fileNumbers;
            if ~isempty(currentFileNum)
                % Exclude the file number that exceed the number of files
                currentFileNum(currentFileNum > size(tempFiles, 1)) = [];
                if isempty(currentFileNum)
                    error(['Unable to find the images with the file numbers you have specified for subject ', num2str(nSub), ' session ', num2str(nSess)]);
                end
            else
                % Use all files
                currentFileNum = (1:size(tempFiles, 1));
            end
            files(nF).name = tempFiles(currentFileNum, :);
            %end
            % End loop over files
        end
        % End loop over sessions
    end
    % End loop over subjects

    % get the count for time points
    diffTimePoints = ipctb_ica_get_countTimePoints(files);

else


    if isempty(inputFile)

        % dialog Title
        dlg_title = [dataTitle, ' data information.'];

        numParameters = 1;

        % define all the input parameters in a structure
        inputText(numParameters).promptString = 'Number of Subjects?';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '1';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'num_subjects';
        inputText(numParameters).enable = 'on';

        numParameters = numParameters + 1;

        inputText(numParameters).promptString = 'Number of Sessions Per Subject?';
        inputText(numParameters).uiType = 'edit';
        inputText(numParameters).answerString = '1';
        inputText(numParameters).dataType = 'numeric';
        inputText(numParameters).tag = 'num_sessions';
        inputText(numParameters).enable = 'on';

        %             numParameters = numParameters + 1;
        %
        %             inputText(numParameters).promptString = 'Select data type.';
        %             inputText(numParameters).uiType = 'popup';
        %             inputText(numParameters).answerString = {'Real', 'Complex'};
        %             inputText(numParameters).dataType = 'string';
        %             inputText(numParameters).tag = 'datatype';
        %             inputText(numParameters).enable = 'on';
        %
        %             numParameters = numParameters + 1;
        %
        %
        %             inputText(numParameters).promptString = 'Select an option to read complex images.';
        %             inputText(numParameters).uiType = 'popup';
        %             inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
        %             inputText(numParameters).dataType = 'string';
        %             inputText(numParameters).tag = 'read_complex_images';
        %             inputText(numParameters).enable = 'off';
        %
        %             numParameters = numParameters + 1;
        %
        %
        %             inputText(numParameters).promptString = 'Select an option to write complex images.';
        %             inputText(numParameters).uiType = 'popup';
        %             inputText(numParameters).answerString = {'Real&Imaginary', 'Magnitude&Phase'};
        %             inputText(numParameters).dataType = 'string';
        %             inputText(numParameters).tag = 'write_complex_images';
        %             inputText(numParameters).enable = 'off';

        numUIControls = length(inputText);

        % Input dialog box (get the necessary numbers)
        answer = ipctb_ica_inputDialog('inputtext', inputText, 'Title', dlg_title);

        if ~isempty(answer)
            % number of subjects
            numOfSub = answer{1};
            % number of sessions
            numOfSess = answer{2};
            getSPMMatrix = 'no'; %answer{3};
            %                 dataType = answer{3};
            %                 % read complex images
            %                 read_complex_images = answer{4};
            %                 % write complex images
            %                 write_complex_images = answer{5};
        else
            error('Number of data sets is not specified');
        end

        fileInfo.fileName = subjectFile;
        fileInfo.format = '';

        %             sesInfo.userInput.read_complex_images = lower(read_complex_images);
        %
        %             sesInfo.userInput.write_complex_images = lower(write_complex_images);


        % for one subject and one session
        if numOfSub*numOfSess == 1
            % select data
            files(1).name = ipctb_ica_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', ...
                'filter', '*.img;*.nii', 'title', 'Select files for subject 1 session 1');
            if isempty(files(1).name)
                error([dataTitle, ' data is not selected for the analysis']);
            else
                % get count for time points
                [diffTimePoints] = ipctb_ica_get_countTimePoints(files);
            end

        else

            % select the data
            [files, diffTimePoints] = ipctb_ica_select_data('title', ['Select ', modalityType, ' data'], 'num_subjects', ...
                numOfSub, 'num_sessions', numOfSess, 'files_specification', 'unequal', 'spm_check', 'no', ...
                'filter_string', '*.img;*.nii', 'type_file_selection', 'multiple', 'fileInfo', fileInfo, 'figure_menu', 'data', ...
                'datatype', dataType, 'complex_file_naming', read_complex_file_naming, ...,
                'read_complex_images', read_complex_images);
        end

    else

        %--number of subjects
        keywd = 'selectedSubjects';
        inputData = ipctb_ica_read_variables(inputFile, keywd, 'vector', 'cell');
        [selectedSubjects] = getfield(inputData, keywd);

        numOfSub = length(selectedSubjects);

        %--number of sessions
        keywd = 'numOfSess';
        inputData = ipctb_ica_read_variables(inputFile, keywd, 'scalar', 'integer');
        numOfSess = getfield(inputData, keywd);
        if numOfSess == 0
            error('Number of sessions cannot be zero');
        end

        disp('Reading data ...');

        % Loop over number of subjects
        for i = 1 : numOfSub
            % Strings corresponding to the selected subjects
            subStr = selectedSubjects{i};

            % one design matrix for each subject
            if strcmpi(spmMatFlag, 'diff_sub_diff_sess')
                % For all subjects search for the suffix designMat
                keywd = [subStr '_designMat'];
                tempInputData = ipctb_ica_read_variables(inputFile, keywd, 'scalar', 'string');
                ipctb_ica_checkVariable(tempInputData, keywd, inputFile);
                designMatrix(i).name = getfield(tempInputData, keywd);
            end

            % loop over sessions
            for j = 1 : numOfSess
                sessStr = ['_s', num2str(j)];
                keywd = [subStr sessStr];
                tempInputData = ipctb_ica_read_variables(inputFile, keywd, 'vector', 'cell');
                ipctb_ica_checkVariable(tempInputData, keywd, inputFile);
                Value_vector = getfield(tempInputData, keywd);
                clear tempInputData;
                sub(i).ses(j).dir = Value_vector{1};
                sub(i).ses(j).imagePattern = Value_vector{2};
                if length(Value_vector) == 3
                    sub(i).ses(j).fileNum = Value_vector{3};
                else
                    sub(i).ses(j).fileNum = [];
                end
            end
            % end for loop over sessions
        end
        % end for loop over subjects

        counter = 0;
        % Initialise input file information
        files = repmat(struct('name', []), 1, numOfSub*numOfSess);
        for j = 1 : numOfSub
            for k = 1 : numOfSess
                counter = counter + 1;
                % get scans and parse
                fileDir = sub(j).ses(k).dir;
                filePattern = sub(j).ses(k).imagePattern;
                fileNum = sub(j).ses(k).fileNum;

                % list Files with the matching pattern
                [fileList] = ipctb_ica_listFiles_inDir(fileDir, filePattern);

                if(isempty(fileList))
                    infoCell{1} = fileDir;
                    infoCell{2} = filePattern;
                    ipctb_ica_error('Could not find any files matching pattern', infoCell);
                end
                % get the full file path for the files
                fileListWithDir = ipctb_ica_fullFile('directory', fileDir, 'files', fileList);
                % checks only the first file
                testFile = deblank(fileListWithDir(1, :));
                if ~exist(testFile, 'file')
                    error(['File ', testFile, ' does not exist']);
                end
                fileListWithDir = ipctb_ica_rename_4d_file(fileListWithDir);
                if ~isempty(fileNum)
                    fileNum(fileNum > size(fileList, 1)) = [];
                    if isempty(fileNum)
                        error(['Unable to find the images with the file numbers you have specified for subject ', num2str(j), ' session ', num2str(k)]);
                    end
                    fileListWithDir = fileListWithDir(fileNum, :);
                end
                files(counter).name = fileListWithDir; % append the file list with directory
            end
        end

        % Count for the time points
        diffTimePoints = ipctb_ica_get_countTimePoints(files);

        disp('Done reading data ...');

    end

end

% Check spm design matrix
ipctb_ica_check_spm_design_matrix(designMatrix, numOfSub, numOfSess, diffTimePoints, spmMatFlag, inputFile);
