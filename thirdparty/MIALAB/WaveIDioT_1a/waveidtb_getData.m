function waveidtb_getData(hObject)

global GUI_GLOBAL_DATA;

% Define the figure window
getData_Handle = figure( 'Visible' , 'off' , 'Position' , [50,200,550,600] ,'NumberTitle','off' , ...
    'MenuBar' , 'none' , 'Toolbar' , 'none', 'Resize' , 'off',...
    'Name','Setup Denoising Analysis','Color','k');
movegui(getData_Handle, 'center');
set(gcf,'Visible','on');
% Title for this window
getData_title = ' Specify Parameters ';
titleColor = 'y';
titleSize = 20;
axes('Parent', getData_Handle, 'position', [0 0 1 1], 'visible', 'off');
text( 0.5, 0.925, getData_title, 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', ...
    'fontsize', titleSize, 'HorizontalAlignment', 'center', 'FontName', 'times');

%% Define the Common Color Properties for this Dialog box
titleColor = 'w';
titleSize = 13;

%% Parameter 1
param1Title = ' 1. Type of Wavelet basis ?  ';
% Text
param1Pos = [0.05 0.65 0.5 0.1];
uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',param1Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',param1Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times');

% Dropdown Menu
select_waveletType_menu = uicontrol('Parent',getData_Handle,'Style','popupmenu',...
    'String',{'Default','db2','coif2','bior1.3'},'Units', 'Normalized',...
    'Position',[.7 .65 .2 .1],'BackgroundColor','k','ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' ,  'Callback',{@select_waveletType_callback});

%% Parameter 2
param2Title = ' 2. Number of decomposition levels ?  ';
% Text
param2Pos = [0.05 0.55 0.6 0.1];
uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',param2Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',param2Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times');

% Dropdown Menu
select_numLevels_menu = uicontrol('Parent',getData_Handle,'Style','popupmenu',...
    'String',{'Default','3','4'},'Units', 'Normalized',...
    'Position',[.7 .55 .2 .1],'BackgroundColor','k','ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' ,  'Callback',{@select_numLevels_callback});

%% Parameter 3
param3Title = ' Have you selected the data ? ';
% Text
param3Pos = [0.05 0.45 0.55 0.1];
select_data_text = uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',param3Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',param3Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times');

% Pushbutton to get the data
select_data_button = uicontrol('Parent',getData_Handle,'Style','pushbutton',...
    'String','Select','Units', 'Normalized',...
    'Position',[.65 .5 .25 0.05],'BackgroundColor',[0 0.5 0.5],'ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' ,  'Callback',{@select_data_callback});

%% Parameter 4 (option 1)
option1Title = ' 3. Select image files using one option: ';
select_option1_text = uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',option1Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',param3Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times' , 'Visible', 'off');

% Pushbutton to get the data
select_option1_menu = uicontrol('Parent',getData_Handle,'Style','popupmenu',...
    'String',{'Default', 'SPM_SELECT GUI','DATA_IN_ONE_FOLDER'},'Units', 'Normalized',...
    'Position',[.7 .45 .2 0.1],'BackgroundColor','k','ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' , 'Callback',@select_option1_callback, 'Visible', 'off');

%% Parameter 5 (option 2)
% Ask PREFIX of Data
option2Title = '4. What is the prefix of data to be searched ?';
option2Pos =  [0.05 0.35 0.55 0.1];
select_option2_text = uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',option2Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',option2Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times' , 'Visible', 'off');
% Enter Text
select_option2_menu = uicontrol('Parent',getData_Handle,'Style','edit',...
    'Units', 'Normalized',...
    'Position',[.7 .4 .20 .05],'BackgroundColor','k','ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' , 'Callback',@select_option2_callback,'Visible', 'off');

%% Parameter 6 (option 3)
% Ask FORMAT of Data
option3Title = '5. What is the data format ?';
option3Pos =  [0.05 0.25 0.55 0.1];
select_option3_text = uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',option3Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',option3Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times' , 'Visible', 'off');
% Enter Text
select_option3_menu = uicontrol('Parent',getData_Handle,'Style','popupmenu',...
    'String',{'Default','NIFTI','ANALYZE'},'Units', 'Normalized',...
    'Position',[.7 .25 .2 .1],'BackgroundColor','k','ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' , 'Callback',@select_option3_callback,'Visible', 'off');

%% Parameter 7 (option 4)
% Ask FORMAT of Data
option3Title = '6. Select root directory for all Subject(s)/Session(s) : ';
option3Pos =  [0.05 0.15 0.55 0.1];
select_option4_text = uicontrol('Parent',getData_Handle,'Style','text','BackgroundColor','k',...
    'String',option3Title,'ForegroundColor',titleColor,'Units','Normalized',...
    'Position',option3Pos,'FontWeight','bold', ...
    'FontSize', titleSize,'HorizontalAlignment','left','FontName','times' , 'Visible', 'off');
% Enter Text
select_option4_menu = uicontrol('Parent',getData_Handle,'Style','pushbutton',...
    'String','Select','Units', 'Normalized',...
    'Position',[.7 .2 .20 0.05],'BackgroundColor',[0 0.5 0.5],'ForegroundColor',...
    titleColor,'FontSize', titleSize,'HorizontalAlignment','Center',...
    'FontName','times' , 'Callback',@select_option4_callback,'Visible', 'off');

%% Done (Pushbutton)

select_done_button = uicontrol('Parent',getData_Handle,'Style','pushbutton',...
    'String','DONE','Units', 'Normalized',...
    'Position',[.7 .05 .2 0.05],'BackgroundColor',[0 0.5 0.5],'ForegroundColor',...
    titleColor,'FontSize', 14 ,'HorizontalAlignment','Center',...
    'FontName','times' ,  'Callback',{@select_done_callback});

set(select_done_button , 'Enable','Off');

%% Define global variables for this window
guidata( getData_Handle , GUI_GLOBAL_DATA );


%% Waitbar
h_wait = waitbar(0,'Please Wait...','CloseRequestFcn',@close_waitbar,'Visible', 'Off');

    function close_waitbar(hObject,eventdata)
        delete(gcbf)
    end

%% Callback Functions

% Callback to Select the Wavelet basis functions (Parameter 1)
    function select_waveletType_callback(source,eventdata)
        TYPE = get(select_waveletType_menu,'Value');
        
        if TYPE==2
            GUI_GLOBAL_DATA.WaveType = 'db2';
            
        elseif TYPE==3;
            GUI_GLOBAL_DATA.WaveType = 'coif2';
            
        elseif TYPE==4;
            GUI_GLOBAL_DATA.WaveType = 'bior1.3';
        else % Default
            GUI_GLOBAL_DATA.WaveType = 'sym2';
        end;
        
        % Update global variable
        guidata( getData_Handle , GUI_GLOBAL_DATA );
        guidata( hObject , GUI_GLOBAL_DATA );
    end


% Callback to Select the Number of wavelet decomp. levels (Parameter 2)
    function select_numLevels_callback(source,eventdata)
        
        NumLev = get(select_numLevels_menu,'Value');
        
        if NumLev==2
            GUI_GLOBAL_DATA.NumLev = 3;
            
        elseif NumLev==3;
            GUI_GLOBAL_DATA.NumLev = 4;
        else % Default
            GUI_GLOBAL_DATA.NumLev = 2;
        end;
        % Update global variable
        guidata( getData_Handle , GUI_GLOBAL_DATA );
        guidata( hObject , GUI_GLOBAL_DATA );
    end


% Callback to get the location of data (using SPM_Select or ) (Parameter 3)
    function select_data_callback(source,eventdata)
        
        % Replace the Select button with yes or no
        set(select_data_button,'Visible','off');
        set(select_data_text,'Visible','off');
        
        % Option 1 (replaces Parameter 3)
        % Create Drop down menu
        set(select_option1_menu,'Visible','on');
        set(select_option1_text,'Visible','on');
        guidata( hObject , GUI_GLOBAL_DATA );
    end

global len;

    function select_option1_callback(source,eventdata)
        
        option1 = get(select_option1_menu,'Value');
        
        if option1 == 2
            
            set(select_option2_text,'Visible','off');
            set(select_option2_menu,'Visible','off');
            set(select_option3_text,'Visible','off');
            set(select_option3_menu,'Visible','off');
            
            [GUI_GLOBAL_DATA.fileNames,GUI_GLOBAL_DATA.dirs] = spm_select(Inf,'image');
            GUI_GLOBAL_DATA.fileNames = mat2cell(GUI_GLOBAL_DATA.fileNames ,...
                ones(size(GUI_GLOBAL_DATA.fileNames,1),1)' , ...
                length(GUI_GLOBAL_DATA.fileNames(1,:)));
            
            len = length(GUI_GLOBAL_DATA.fileNames{1});
            GUI_GLOBAL_DATA.fileFormat = GUI_GLOBAL_DATA.fileNames{1}((len-4):(len-2));
            guidata( getData_Handle , GUI_GLOBAL_DATA );
            set(select_done_button , 'Enable','On');
            GUI_GLOBAL_DATA.dataSelectOption = option1;
        elseif option1 == 3;
            set(select_option2_text,'Visible','on');
            set(select_option2_menu,'Visible','on');
            GUI_GLOBAL_DATA.dataSelectOption = option1;
        else % Default;
            [GUI_GLOBAL_DATA.fileNames,GUI_GLOBAL_DATA.dirs] = spm_select(Inf,'image');
            set(select_done_button , 'Enable','On');
        end;

        guidata( getData_Handle , GUI_GLOBAL_DATA );
        guidata( hObject , GUI_GLOBAL_DATA );
    end


    function select_option2_callback(source,eventdata)
        %             waitfor(select_option2_menu , 'String');
        GUI_GLOBAL_DATA.origDataPrefix = get(select_option2_menu,'String');
        guidata( getData_Handle , GUI_GLOBAL_DATA );
        guidata( hObject , GUI_GLOBAL_DATA );
        set(select_option3_menu,'Visible','on');
        set(select_option3_text,'Visible','on');
    end


    function select_option3_callback(source,eventdata)
        
        
        if isempty(GUI_GLOBAL_DATA.origDataPrefix)
            msgbox('Please enter the prefix of your image files. If they have different prefix, then please SPM_SELECT GUI option to manually select the data separately','INCOMPLETE INFORMATION!!!','error');
            guidata( getData_Handle , GUI_GLOBAL_DATA );
%             select_option2_callback;
        end;

        % Get the File format Value from the GUI (popupmenu)
        GUI_GLOBAL_DATA.fileFormat = get(select_option3_menu,'Value');
        if GUI_GLOBAL_DATA.fileFormat == 2
           
            GUI_GLOBAL_DATA.fileFormat = 'nii';
           
       elseif GUI_GLOBAL_DATA.fileFormat == 3

           GUI_GLOBAL_DATA.fileFormat = 'img';
           
        else
            
            GUI_GLOBAL_DATA.fileFormat = 'nii';
            
        end;
        
        set(select_option4_menu,'Visible','on');
        set(select_option4_text,'Visible','on');
        guidata( getData_Handle , GUI_GLOBAL_DATA );
        guidata( hObject , GUI_GLOBAL_DATA );
    end

    function select_option4_callback(source,eventdata)
    
        % Acquire the data using option1 = 3 (Check Again)
        if GUI_GLOBAL_DATA.dataSelectOption == 3;
            prefix =  GUI_GLOBAL_DATA.origDataPrefix;
            fileFormat =  GUI_GLOBAL_DATA.fileFormat;
            
            if isunix
                 set(h_wait,'Visible','on');
                [ dirname ] = uigetdir('/homes/','Select data directory (root dir for all subjects)');
            else
                 set(h_wait,'Visible','on');
                [ dirname ] = uigetdir('C:\','Select data directory (root dir for all subjects)');
            end;
            if dirname
                set(h_wait,'Visible','on');
                waitbar(1/2,h_wait,'Please wait... Reading image paths from directory..');
                searchVariable = ['^',char(prefix),'.*\.',fileFormat,'$'];
                [GUI_GLOBAL_DATA.fileNames,GUI_GLOBAL_DATA.dirs] = spm_select('ExtFPListRec', dirname , searchVariable , GUI_GLOBAL_DATA.timePoints);
                
                waitbar(1,h_wait,'Done');
                set(h_wait,'Visible','off');
                set(select_done_button , 'Enable','On');
            else
                msgbox('No data directory selected! Please reselect option 5.','error');
            end;
        else
            set(select_option1_menu,'Visible','off');
            set(select_option1_text,'Visible','off');            
            set(select_option2_menu,'Visible','off');
            set(select_option2_text,'Visible','off');            
            set(select_option3_menu,'Visible','off');
            set(select_option3_text,'Visible','off');
            set(select_option4_menu,'Visible','off');
            set(select_option4_text,'Visible','off');
        end;
        guidata( getData_Handle , GUI_GLOBAL_DATA );
        guidata( hObject , GUI_GLOBAL_DATA );
    end

    function select_done_callback(source,eventdata)
        
        % Save MAT file to memory
        if strcmpi( GUI_GLOBAL_DATA.outDirFlag , 'UserDefined' );
            GUI_GLOBAL_DATA.outDir = GUI_GLOBAL_DATA.outDir;
        else strcmpi( GUI_GLOBAL_DATA.outDirFlag , 'Default' );
            if isunix
                GUI_GLOBAL_DATA.outDir = '/homes/';
            else
                GUI_GLOBAL_DATA.outDir = pwd;
            end;
        end;
        GUI_GLOBAL_DATA.ParameterMATfileName = fullfile(GUI_GLOBAL_DATA.outDir, ['waveidtb_ParameterInfo_',datestr(now,'ddmmmyyyy_HHMMSS')]);
        save(GUI_GLOBAL_DATA.ParameterMATfileName, 'GUI_GLOBAL_DATA');
        guidata( getData_Handle , GUI_GLOBAL_DATA );
        delete(gcbf);
        fprintf('Denoising Paremeters were stored to memory.\n');
        disp(['Filename stored (with path): ',GUI_GLOBAL_DATA.ParameterMATfileName,'.\n']);
        fprintf('Continue to "Run Denoising".\n');
        guidata( hObject , GUI_GLOBAL_DATA );
    end

guidata(hObject,GUI_GLOBAL_DATA);

end