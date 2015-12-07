function outmat = zpad(inmat,dim,lvl)
%
% FUNCTION:     zpad    -- Pads zeros to input 3-D matrix in x,y or z 
%                          direction.
%                            
% USAGE:        outmat = zpad(inmat,dim,lvl)
%
% Inputs:       inmat   -- 3-D matrix input for zero padding.
%               dim     -- Dimension for circularizing. Input must be a
%                          string. Valid dimensions: 'x','y','z'.
%               lvl     -- Pad level - how many columns to zeros to pad.
%
% Output:       outmat  -- Zero-padded matrix.
%
% EXAMPLE:      output = zpad(3Dmat,'x',2)
%               
%               3Dmat (:,:,1)   = [a b . . . . . c d]    
%                                 [a b . . . . . c d]
%                                 [a b . . . . . c d]
%                                 [a b . . . . . c d]    
%
%               output (:,:,1)  = [0 0 a b . . . . . c d 0 0]
%                                 [0 0 a b . . . . . c d 0 0]
%                                 [0 0 a b . . . . . c d 0 0]
%                                 [0 0 a b . . . . . c d 0 0]
%
% AUTHOR:       Ameera X Patel
% CREATED:      29-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     3
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: zpad.m 3 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

valid_dims={'x','y','z'};

err=struct();
err(1).inp=chkninput(nargin,[3,3],nargout,[0,1],fname);

intype=class(inmat);
if ~(strcmpi(intype,'double') || strcmpi(intype,'logical'))
    err(2).inp=1;
else
    err(2).inp=0;
end

err(3).inp=chkintype(dim,'char',fname,valid_dims);
err(4).inp=chkintype(lvl,'numeric',fname);

if sum(cat(1,err.inp))>=1
    outmat=[];
    return
end

%% compute dimensions and calculate pad level

[R,C,N]=size(inmat);

lvl=round(lvl);

%% zeropad in x,y or z dimension

if      strcmpi(dim,'x')
        outmat=zeros(R,(C+(lvl*2)),N);
        outmat(:,(1+lvl:C+lvl),:)=inmat;
        
elseif  strcmpi(dim,'y')
        outmat=zeros((R+(lvl*2)),C,N);
        outmat((1+lvl:R+lvl),:,:)=inmat;
        
elseif  strcmpi(dim,'z')
        outmat=zeros(R,C,(N+(lvl*2)));
        outmat(:,:,(1+lvl:N+lvl))=inmat;
end

end