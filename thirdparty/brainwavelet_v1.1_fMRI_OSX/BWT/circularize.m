function outmat = circularize(inmat,dim,lvl)
%
% FUNCTION:     circularize -- Circularizes 3-D matrix in x,y or z
%                              direction.
%                            
% USAGE:        outmat = circularize(inmat,dim,lvl)
%
% Inputs:       inmat       -- 3-D matrix input for circularizing.
%               dim         -- Dimension (valid dimensions: x','y','z') for 
%                              circularizing. Input must be a string.
%               lvl         -- Pad level - how many columns to circularize.
%
% Output:       outmat      -- Circularized matrix.
%
% EXAMPLE:      output = circularize(3Dmat,'x',2)
%               
%               3Dmat (:,:,1)   = [a b . . . . . c d]    
%                                 [a b . . . . . c d]
%                                 [a b . . . . . c d]
%                                 [a b . . . . . c d]    
%
%               output (:,:,1)  = [c d a b . . . . . c d a b]
%                                 [c d a b . . . . . c d a b]
%                                 [c d a b . . . . . c d a b]
%                                 [c d a b . . . . . c d a b]
%
% AUTHOR:       Ameera X Patel
% CREATED:      29-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: circularize.m 4 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return;
end

err=struct();
err(1).inp=chkninput(nargin,[3,3],nargout,[0,1],fname);

intype=class(inmat);
if ~(strcmpi(intype,'double') || strcmpi(intype,'logical'))
    err(2).inp=1;
else
    err(2).inp=0;
end

err(3).inp=chkintype(dim,'char',fname,{'x','y','z'});
err(4).inp=chkintype(lvl,'numeric',fname);

if sum(cat(1,err.inp))>=1;
    outmat=[]; return;
end

clear fname err

%% compute dimensions and calculate circ pad

[R,C,N]=size(inmat);

lvl=round(lvl);

%% circularize in x,y or z dimension

if      strcmpi(dim,'x')
        outmat=zeros(R,(C+(lvl*2)),N);
        outmat(:,(1:lvl),:)=inmat(:,((C-(lvl-1)):C),:);
        outmat(:,(1+lvl:C+lvl),:)=inmat;
        outmat(:,((C+lvl+1):(C+(lvl*2))),:)=inmat(:,(1:lvl),:);
        
elseif  strcmpi(dim,'y')
        outmat=zeros((R+(lvl*2)),C,N);
        outmat((1:lvl),:,:)=inmat(((R-(lvl-1)):R),:,:);
        outmat((1+lvl:R+lvl),:,:)=inmat;
        outmat(((R+lvl+1):(R+(lvl*2))),:,:)=inmat((1:lvl),:,:);

elseif  strcmpi(dim,'z')
        outmat=zeros(R,C,(N+(lvl*2)));
        outmat(:,:,(1:lvl))=inmat(:,:,((N-(lvl-1)):N));
        outmat(:,:,(1+lvl:N+lvl))=inmat;
        outmat(:,:,((N+lvl+1):(N+(lvl*2))))=inmat(:,:,(1:lvl));
end

end