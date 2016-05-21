function vox=mm2vox(varargin)
% function vox=mm2vox(mm, V)
%
% converts mm-coordinates in voxel-coordinates/indices
% using an existing transformation matrix 
% vox = [x, y, z] voxel coordinates
% mm  = [x y z] mm coordinates
% V = structure containing Image information (e.g., use the MNI template) -
% (see spm_vol)

% uses spm_get_space.m by John Ashburner 
% 

if nargin<1
    mfile_showhelp; 
    return;
else
    mm= varargin{1};
    d = find(size(mm) == 3);
    if isempty(d),
        error('wrong usage: vox = mm2vox(mm-coord, [fname])');
    elseif d==1 && length(d)==1,
        mm=mm'; %arrange coordinate triples in rows        
    end
    
    if nargin<2 
        if exist('spm_select','file')
            fname=spm_select(1,'image','select file defining space');
        else
            fname=spm_get(1,'*.img','select file defining space');
        end    
        V=spm_vol(fname);
    else
        if isstruct(varargin{2}), 
            V=varargin{2};
        else
            V(1).fname=varargin{2};
        end    
    end
end

v2m=spm_get_space(V(1).fname);
m2v=inv(v2m);

for i=1:size(mm,1)
    vox(i,1:3)=mm(i,:)*m2v(1:3,1:3) + m2v(1:3,4)';
end    
