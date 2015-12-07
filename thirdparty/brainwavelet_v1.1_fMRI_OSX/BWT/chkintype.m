function x = chkintype(var,type,rootfname,strmatch)
%
% FUNCTION:     chkintype.m -- Error function that checks type of input 
%                              variable (and string matching options).
%                
% USAGE:        x = chkintype(var,type,rootfname,strmatch)
%
% Inputs:       var         -- Number of actual inputs.
%               type        -- Limits of accepted number of inputs.
%               rootfname   -- Filename that called this function.
%               strmatch    -- Accepted values. Input must be a string.
%
% Output:       x           -- Error counter.
%
% EXAMPLE:      errors = chkinput(variable,'char',mfilename,...
%                        [cat,dog,bat,rat])
%
% AUTHOR:       Ameera X Patel
% CREATED:      23-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     5
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: chkintype.m 5 03-02-2014 BWTv1.1 axpatel

                    
%% check number of input arguments to this function

currfname=mfilename;
if nargin<1
    help(currfname); return;
end
x=0;

%narginchk(3,4)
%nargoutchk(0,1)

%% check inputs are of the right format

if ~ischar(rootfname)
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Error using %s: function name must be a string \n\n',...
        currfname);
    x=x+1;
    return
end

if ~ischar(type)
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Error using %s in %s: input type must be a string \n\n',...
        currfname,rootfname);
    x=x+1;
    help chkintype
    return
end

%% validate fundamental types

validtypes={'single','double','numeric','float',...
            'int8','int16','int32','int64',...
            'uint8','uint16','uint32','uint64'...
            'logical','char','struct','cell','function_handle'};
        

if sum(strcmpi(type,validtypes))~=1
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Error using %s in %s: invalid variable type \n\n',...
        currfname,rootfname,typestr);
    x=x+1;
    return
end

%% check variable type is correct

if isa(var,type)==0
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Error using %s in %s: input variable must\nbe of type %s\n\n',...
        currfname,rootfname,type);
    x=x+1;
end

%% match variable string to valid type

ok_listtype={'single','double','numeric','char','struct','cell'};
ok_vartype={'single','double','numeric','char'};

% convert strmatch type is ok and convert to cell array for matching

if exist('strmatch','var')
    if sum(strcmpi(class(strmatch),ok_listtype))==0
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','Error using %s in %s: list for matching\nmust be cell/numeric array or struct\n\n',...
            currfname,rootfname);
        x=x+1;
        return
    end
    if isnumeric(strmatch)
        strmatch=arrayfun(@num2str,strmatch,'unif',0);
    end
    if ischar(strmatch)
        strmatch=cellstr(strmatch);
    end
    if isstruct(strmatch)
        strmatch=cellfun(@num2str,struct2cell(strmatch),'unif',0);
    end
    
    if ~iscell(strmatch)
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','Error using %s in %s: list for matching\nmust be cell/numeric array or struct\n\n',...
            currfname,rootfname);
        x=x+1;
        return
    end
    
    if ~ischar(var)
        if sum(strcmpi(class(var),ok_vartype))==1
            var=num2str(var);
        else
            cprintf('_err','*** BrainWavelet Error ***\n')
            cprintf('err','Error using %s in %s: cannot parse\ninput variable type - type must be string or number\n\n',...
                currfname,rootfname);
            x=x+1; 
        end
    end
    
    % compare var and strmatch
    
    if sum(strcmpi(var,strmatch))==0;
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','Error using %s in %s: value of input\nvariable is not allowed\n\n',...
                currfname,rootfname);
        x=x+1; 
        return
    end
end

end