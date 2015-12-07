function [] = header
%                  
% FUNCTION:     header -- Wavelet despiking header function.
%                
% USAGE:        header()
%
% AUTHOR:       Ameera X Patel
% CREATED:      10-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     8
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: header.m 8 27-03-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;

if chkninput(nargin,[0,0],nargout,[0,0],fname)>=1;
    return
end

%% header info
cprintf([0.1,0.1,0.7],'\n============================================\n')
cprintf([0.1,0.1,0.7],'         BrainWavelet Toolbox v1.1')
cprintf([0.1,0.1,0.7],'\n============================================\n')
%cprintf([1,0.4,0],'\n------------------------------------------------\n')
%cprintf([1,0.4,0],'************* ')
%cprintf([0.1,0.1,0.7],'BrainWavelet Toolbox ')
%cprintf([1,0.4,0],'*************\n')
%cprintf([1,0.4,0], '------------------------------------------------\n')

cprintf([0,0,0],'\nAuthor: ')
cprintf([0.1,0.1,0.7],'Ameera Patel, 2014\n')

%cprintf([0,0,0], '\nPlease cite:  ')
%cprintf([0.1,0.1,0.7],'Patel et al., A wavelet method for \n')
%cprintf([0.1,0.1,0.7],'modeling and denoising motion artefacts from \n')
%cprintf([0.1,0.1,0.7],'resting-state fMRI time series.\n')
%cprintf([0.1,0.1,0.7],'Neuroimage [in peer review].\n\n')

%cprintf([1,0.4,0],'------------------------------------------------\n\n')
fprintf('\nWavelet Despiking time series ...\n\n')

end