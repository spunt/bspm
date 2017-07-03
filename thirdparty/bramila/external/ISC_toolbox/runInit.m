% When using scripts, keep the following directories up to date:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis destination directory. This must be same as
% specified in the parameter file.
analysisDestinationPath = '/share/sig_pic/pajula/ISCtesting/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory of the templates:
templatePath = '/share/sig_pic/pajula/ISCToolbox/templates/';
% Directory of the codes for running the ISC analysis:
ISCtoolboxPath = '/share/sig_pic/pajula/ISCofficial/isc-toolbox/';
% Directory of the codes to process nifti -files:
niftiPath = '/share/sig_pic/pajula/ISCofficial/isc-toolbox/niftitools/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not edit this part:
if exist(templatePath) ~= 7
    warning(['MNI templates not found in ' templatePath '.'])
end
if exist(ISCtoolboxPath) ~= 7
    error(['ISC toolbox not found in ' ISCtoolboxPath '.'])
end
if exist(templatePath) ~= 7
    warning(['Destination directory ' analysisDestinationPath ' does not exist.'])
end
curDir = cd;
cd(ISCtoolboxPath) 
setISCToolboxPath;
cd(curDir)

% set path:
disp(' ');disp(' ')
P=path;
P=path(P,ISCtoolboxPath);
P=path(P,analysisDestinationPath);
P=path(P,niftiPath);
P=path(P,templatePath);

disp(' ');disp(' ')
% load parameters from the destination directory. 
% Then you can run analysis by typing: runAnalysis(Params).
try
  load([analysisDestinationPath 'Tag'])
  load([analysisDestinationPath Tag '.mat'])
  [Params,flag] = paramsValidation(Params);
  if ~flag
    disp('Current parameters loaded to workspace.');disp(' ')
    disp('If you fill in parameters manually, validate parameters by typing: "Params=paramsValidation(Params);"');disp(' ')
    disp('Type "ISCanalysis" to use GUI to set parameters.');disp(' ')    
  end            
catch
    warning(['No parameter struct found in ' analysisDestinationPath]);disp(' ')
    disp('Default parameters loaded to workspace.');
    Params = initParams(templatePath,analysisDestinationPath);
    [Params,flag] = paramsValidation(Params);    
    if ~flag
        disp('If you fill in parameters manually, Type "Params.PublicParams" to see them.');disp(' ')
        disp('Validate parameters by typing: "Params=paramsValidation(Params);"');disp(' ')
        disp('Type "ISCanalysis" to use GUI to set parameters!')
    end
end

disp(' ')

if ~flag
    error('Cannot proceed to analysis!')
end
