function [xMNI,yMNI,zMNI] = bramila_space2MNI(x,y,z)
% BRAMILA_SPACE2MNI - Converts matrix indeces into MNI 2mm space. At the moment this only works for
% the FSL pipeline i.e. volumes of size [91 109 91]
% - Usage:
%	[xMNI, yMNI, zMNI] = bramila_space2MNI(x,y,z);

% EG 2014-01-14

    xMNI = 2*((x-1)-45);
    yMNI = 216*(y-1)/108 - 126;
    zMNI = 180*(z-1)/90 - 72;
end   
