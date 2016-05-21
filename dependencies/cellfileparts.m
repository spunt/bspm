function [pth, fname, fext] = cellfileparts(incell)
% CELLFILEPARTS Fileparts for cell arrays
%
%  USAGE: [pth, fname, fext] = cellfileparts(incell)
%
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-06-26
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
if nargin==0, mfile_showhelp; return; end
[pth, fname, fext] = cellfun(@fileparts, incell, 'Unif', false); 