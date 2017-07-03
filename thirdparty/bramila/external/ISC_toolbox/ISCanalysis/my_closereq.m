% my_closereq
% User-defined close request function
% to display a question dialog box

selection = questdlg('Exit LSA tool?',...
    'Close Request Function',...
    'Yes','No','Yes');
switch selection,
    case 'Yes',
        delete(gcf)
    case 'No'
        return
end