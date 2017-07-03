function [x,y,z] = bramila_MNI2space(xMNI,yMNI,zMNI)
% BRAMILA_MNI2SPACE - Converts 2mm MNI coordinates into matrix indeces. At the moment this only works for
% the FSL pipeline i.e. volumes of size [91 109 91]
% - Usage:
%	[x,y,z] = bramila_MNI2space(xMNI,yMNI,zMNI);

% EG 2014-01-14
% Notes: add support for more millimeters sizes 


    x = xMNI/2 + 46;
    y = 108*(yMNI + 126)/216 +1;
    z=90*(zMNI + 72)/180 +1;
end