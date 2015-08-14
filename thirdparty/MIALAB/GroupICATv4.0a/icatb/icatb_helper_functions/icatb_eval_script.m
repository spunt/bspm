function inputData = icatb_eval_script(file_name)
%% Evaluate script and store the variables in a data structure
%
% Inputs:
% file_name - File name
%
% Outputs:
% inputData - Input data structure
%

oldDir = pwd;

%% Do file parts
[pathstr, fName, extn] = fileparts(file_name);
if isempty(pathstr)
    pathstr = pwd;
end

cd(pathstr);

%% Evaluate file
eval(fName);

cd(oldDir);

vars = whos;

if (length(vars) == 1)
    error(['No variables found in file ', file_name]);
end

% Generate inputData
inputData = struct;
for n = 1:length(vars)
    inputData = setfield(inputData, vars(n).name, eval(vars(n).name));
end
