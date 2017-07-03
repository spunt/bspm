function Params = testParamsMat(Params)
% function tests if all fields exits in the Params struct
% This function maintains the backward compatibility of ISCToolbox
% In the case of older analysis the setup file Params struct does not 
% contain all requested fields of the newer version
%
% Changes in Params-struct over versions:
%
% Version 1.2 (2012), added field for gridComputing and two removal fields
% Version 1.1 (2011), original struct 
% Version 1.0 (2010), original struct 

%Priv = Params.PrivateParams;
Pub = Params.PublicParams;

disp('Testing the compatibility of the loaded Parameters...')

if(~isfield(Pub,'disableGrid'))
    disp('Field disableGrid not found: adding the field...')
    Params.PublicParams.disableGrid = ispc; %If Windows, the grid is disabled
end

if(~isfield(Pub,'removeMemmaps'))
    disp('Field removeMemmaps not found: adding the field...')
    Params.PublicParams.removeMemmaps = false;
end

if(~isfield(Pub,'removeFiltermaps'))
        disp('Field removeFiltermaps not found: adding the field...')
    Params.PublicParams.removeFiltermaps = false;
end
