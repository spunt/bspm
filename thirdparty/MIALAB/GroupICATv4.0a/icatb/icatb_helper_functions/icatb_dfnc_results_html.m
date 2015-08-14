function icatb_dfnc_results_html(dfncInfo)
%% Print dfnc results to HTML file
%

helpStr = 'Creating HTML file. This will involve writing jpeg files to the disk. Please wait ...';
helpH = helpdlg(helpStr);
disp(helpStr);

html_file = fullfile(dfncInfo.outputDir, 'html', [dfncInfo.prefix, '_results.html']);
start_string = '<html><head><title> DFNC Results </title></head>';


results = icatb_dfnc_results(dfncInfo, 'fnc oscillations');

titleStr = 'FNC Oscillations';

results_string1 = ['<h2 align = "center">', titleStr, '</h2>'];
for nR = 1:length(results);
    results_string1 = [results_string1, '<p align = "center">', results(nR).text, '</p>'];
    results_string1 = [results_string1, '<p align = "center"> <img src = "', results(nR).file, '" > </img>'];
end

end_string =  '</html>';

results = icatb_dfnc_results(dfncInfo, 'clusters');

titleStr = 'Cluster stats';

results_string2 = ['<<h2 align = "center">', titleStr, '</h2>'];
for nR = 1:length(results);
    results_string2 = [results_string2, '<p align = "center">', results(nR).text, '</p>'];
    results_string2 = [results_string2, '<p align = "center"> <img src = "', results(nR).file, '" > </img>'];
end


results_string = [start_string,  results_string1, results_string2, end_string];

dlmwrite(html_file, results_string, '');

icatb_openHTMLHelpFile(html_file);

try
    delete(helpH);
catch
end

fprintf('Done\n');
