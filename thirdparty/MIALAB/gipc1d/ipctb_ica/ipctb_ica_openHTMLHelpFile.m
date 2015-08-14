function ipctb_ica_openHTMLHelpFile(fileName)
% opens the specified HTML help file in the browser

setupICAFilePath = which(fileName); % get the full file name
[pathstr, fileName, extn] = fileparts(setupICAFilePath); % get the path for the html file
fullFileName = fullfile(pathstr, [fileName, extn]); % form the full file path

s = web(fullFileName, '-browser'); % open HTML file in browser

if s ~= 0
    s = web(fullFileName, '-new');
end