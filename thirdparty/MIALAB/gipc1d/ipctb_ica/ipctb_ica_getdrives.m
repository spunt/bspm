function ret = ipctb_ica_getdrives(varargin)
%GETDRIVES  Get the drive letters of all mounted filesystems on the computer.
%   F = GETDRIVES returns the roots of all mounted filesystems on the computer
%   as a cell array of char arrays.  For UNIX this list consists solely of the
%   root directory, /.  For Microsoft Windows, it is a list of the names of all
%   one-letter mounted drive letters.
%   F = GETDRIVES('-nofloppy') does not scan for the a: or b: drives on
%   Windows, which usually results in annoying grinding sounds coming from
%   the floppy drive.
%   F = GETDRIVES('-twoletter') scans for both one- AND two-letter drive
%   names on Windows.  While two-letter drive names are technically supported,
%   their presence is in fact rare, and slows the search considerably.
%
%   Note that only the names of MOUNTED volumes are returned.  In particular,
%   removable media drives that do not have the media inserted (such as an
%   empty CD-ROM drive) are not returned.
%
%   See also EXIST, COMPUTER, UIGETFILE, UIPUTFILE.

%   Copyright 2001 Bob Gilmore.
%   Email bug reports and comments to bgilmore@mathworks.com

twoletter = logical(0);
nofloppy = logical(0);

% Interpret optional arguments
for i = 1:nargin
    if strcmp('-nofloppy', varargin{i}), nofloppy = logical(1); end
    if strcmp('-twoletter', varargin{i}), twoletter = logical(1); end
end

if ~ispc
    %ret = {'/'};
    ret = {filesep};
else
    % Initialize return cell array
    ret = {};
    
    % Handle -nofloppy flag, or lack thereof.
    startletter = 'a';
    if nofloppy
        startletter = 'c';
    end
    
    for i = double(startletter):double('z')
        if exist([i ':\']) == 7
            ret{end+1} = [i ':\'];
        end
    end
    
    % Handle two-letter case.  The outer loop of this routine could have been
    % combined with the one above, but was left separate for clarity's sake.
    if twoletter
        for i = double('a'):double('z')
            for j = double('a'):double('z')
                if exist([i j ':\']) == 7
                    ret{end+1} = [i j ':\'];
                end
            end
        end
    end
end