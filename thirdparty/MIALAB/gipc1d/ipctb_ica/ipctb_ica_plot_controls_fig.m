function [InputHandle] = ipctb_ica_plot_controls_fig(inputText, figureTag, handle_visibility, okTag, cancelTag, ...
    tagsTobePlotted, InputHandle)
% Plot on a figure the user interface controls with Done and Cancel buttons
% Function callbacks can be defined outside of the function
% 
% Inputs:
% 1. inputText - Structure containing the necessary parameters for plotting
% user interface controls
% 2. figureTag - Figure tag
% 3. handle_visibility - Options are 'on' or 'off'
% 4. okTag - Tag for ok button
% 5. cancelTag - Tag for cancel button
% 6. tagsTobePlotted - User interface controls that are to be plotted.
% 7. InputHandle - Handle for the input figure. If the handle doesn't exist
% new figure will be plotted.
%
% Output:
% 1. InputHandle - Handle for the input figure. If the handle doesn't exist
% new figure will be plotted.


if ~exist('handle_visibility', 'var')
    handle_visibility = 'on';
end

if ~exist('figureTag', 'var')
    figureTag = 'input_dialog';
end

if ~exist('okTag', 'var')
    okTag = 'done';
end

if ~exist('cancelTag', 'var')
    cancelTag = 'cancel';
end

if ~exist('tagsTobePlotted', 'var')
    tagsTobePlotted = cellstr(str2mat(inputText.tag));
end

ipctb_ica_defaults;
global UI_FS;

promptPrefix = 'prompt';

if ~exist('InputHandle', 'var')

    % Setup figure for GUI
    [InputHandle] = ipctb_ica_getGraphics(figureTag, 'normal', figureTag, handle_visibility);

    % Set no menu bar for the figure
    set(InputHandle, 'Menubar', 'none');

end


if ~isfield(inputText, 'value')
    for ii = 1:length(inputText)
        inputText(ii).value = 1;
    end
end

% plot all the user interface controls and their options to the right
% plot a help button to the right

% offsets for x and y
xOffset = 0.02; yOffset = 0.05;

% UI control width and height
uicontrol_Width = 0.3; uicontrol_Height = 0.05;

% number of UIcontrols excluding the push buttons
numUIcontrols = length(tagsTobePlotted);

% Define the push button positions

% Ok Button push button position
okPos(3) = 0.15; okPos(4) = 0.05;
okPos = [0.75 - 0.5*okPos(3) 0.94 - 0.5*okPos(4) okPos(3) okPos(4)];

% Cancel button position
cancelPos = [0.25 - 0.5*okPos(3) okPos(2) okPos(3) okPos(4)];

% set the yPos
yPos = okPos(2) - 1.5*yOffset;

% Position of prompt string
promptPos = [0.02 yPos  0.52 0.1];

% Position of answer string
answerPos = promptPos;
answerPos(1) = promptPos(1) + promptPos(3) + xOffset; answerPos(3) = uicontrol_Width; answerPos(4) = uicontrol_Height;

%%%%%%%%%%% plot all the uicontrols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countTag = 0;
for ii = 1:numUIcontrols
    matchIndex = strmatch(tagsTobePlotted{ii}, str2mat(inputText.tag), 'exact');

    if isfield(inputText, 'uiPos')
        answerPos(3) = inputText(matchIndex).uiPos(1);
        answerPos(4) = inputText(matchIndex).uiPos(2);
    end

    % increment count
    countTag = countTag + 1;
    tempH = ipctb_ica_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
        promptPos, 'string', inputText(matchIndex).promptString, 'tag', [promptPrefix, inputText(matchIndex).tag], 'panel_tag', ...
        ['panel_text', inputText(matchIndex).tag], 'fontsize', UI_FS - 1);

    % get the last handle of tempH
    promptHandle = tempH(end);

    if ~iscell(inputText(matchIndex).promptString)
        newTextString = {inputText(matchIndex).promptString};
    else
        newTextString = inputText(matchIndex).promptString;
    end

    % wrap the prompt string and get the new position
    [newTextString, newPos] = textwrap(promptHandle, newTextString);

    % select the same height as the text wrap
    promptPos(4) = newPos(4);

    promptPos(2) = newPos(2) - 0.5*promptPos(4);
    answerPos(2) = newPos(2) - 0.5*answerPos(4);


    % include panel only for scrolling
    storeTag{countTag} = get(tempH(1), 'tag');

    % store all the initial positions
    initialYPositions(countTag) = promptPos(2);

    set(promptHandle, 'string', newTextString);

    set(tempH(1), 'position', promptPos);

    clear tempH;

    % store all the initial positions
    countTag = countTag + 1;

    tempH = ipctb_ica_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', inputText(matchIndex).uiType, 'position', ...
        answerPos, 'string', inputText(matchIndex).answerString, 'tag', inputText(matchIndex).tag, ...
        'enable', inputText(matchIndex).enable, 'panel_tag', ['panel', inputText(matchIndex).tag], 'value', inputText(matchIndex).value, ...
        'fontsize', UI_FS - 1, 'HorizontalAlignment', 'center');

    storeTag{countTag} = get(tempH(1), 'tag');

    initialYPositions(countTag) = answerPos(2);

    promptPos(2) = promptPos(2) - 0.5*promptPos(4) - yOffset;
    answerPos(2) = answerPos(2) - 0.5*answerPos(4) - yOffset;

end
%%%%%%%%%%% End for plotting all the uicontrols %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Plot push buttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define push buttons
cancelHandle = ipctb_ica_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    cancelPos, 'string', cancelTag, 'Tag', cancelTag);

% define push buttons
okHandle = ipctb_ica_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    okPos, 'string', okTag, 'Tag', okTag);

checkPosition = promptPos(2);

% define slider data
sliderData.tag = storeTag;
sliderData.maxHeight = okPos(2); % - 0.5*okPos(4);
sliderData.minHeight = checkPosition - yOffset;
sliderData.initialYPositions = initialYPositions;

% Plot slider if necessary
if initialYPositions(end) < 0.01

    sliderPos = [0.965 0 0.04 1];

    % Plot slider
    maxVal = 0;
    minVal = -abs(sliderData.minHeight);
    % control the slider step size
    sliderStep = [0.05 0.2];

    sliderH = ipctb_ica_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'slider', 'position', ...
        sliderPos, 'min', minVal, 'max', maxVal, 'value', maxVal, 'callback', {@verSliderCallback, InputHandle}, ...
        'userdata', sliderData, 'sliderstep', sliderStep);

end

%%%%%%%%%%%%%%% Define callbacks for the controls %%%%%%%%%%%%%

function verSliderCallback(handleObj, evd, handles)
% vertical slider callback

% execute the slider callback
ipctb_ica_verSliderCallback(handleObj, handles);