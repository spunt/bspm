function outm = imageToMesh( in,varargin )
%IMAGETOMESH Converts an image (2D) into a mesh (2D) with the image values
%associated to the vertices as an attribute, and a boolean saying weather a
%point is border or not also as second attribute
%
% Optional arguments:
%   'subdivide' adds vertices which improve the mesh shape (divides each square in four triangles)


subdivide=false;
interpolation='NN';

for i=1:numel(varargin)
    if strcmp(varargin{i},'subdivide')
        subdivide=true;
    elseif strcmp(varargin{i},'interpolation')
        interpolation=varargin{i+1};
        i=i+1;
    end
end

nptsIn= prod(in.size);
points_in= in.GetPosition(1:nptsIn)';
allpoints =points_in ;

if subdivide
    
    [ixe,iye] = ndgrid(1:in.size(1)-1,1:in.size(2)-1);
    indices1De=sub2ind(in.size',ixe(:),iye(:))';
    extrapoints = (in.GetPosition(indices1De)+in.spacing*ones(1,numel(indices1De))/2)';
    
    % create topology
    indices1De=sub2ind((in.size-1)',ixe(:),iye(:))';
    neighs=[0 0; 1 0; 1 1 ;0 1; 0 0];
    topology =[];
    for i=1:4
        
        % vertex 2
        [ix,iy] = ndgrid(1+neighs(i,1):in.size(1)-1+neighs(i,1),1+neighs(i,2):in.size(2)-1+neighs(i,2));
        indices1D1=sub2ind(in.size',ix(:),iy(:))';
        % vertex 3
        [ix,iy] = ndgrid(1+neighs(i+1,1):in.size(1)-1+neighs(i+1,1),1+neighs(i+1,2):in.size(2)-1+neighs(i+1,2));
        indices1D2=sub2ind(in.size',ix(:),iy(:))';
        
        topology = [topology ; indices1De(:)+nptsIn indices1D1(:) indices1D2(:) ];
        
        %     [ex,ey] [ix, iy] [ix+1, iy]
        %     [ex,ey] [ix+1, iy] [ix+1, iy+1]
        %     [ex,ey] [ix+1, iy+1] [ix, iy+1]
        %     [ex,ey] [ix, iy+1] [ix, iy]
    end
    
    borders = [1:in.size(1)       1:in.size(1) ones(1,in.size(2))  in.size(1)*ones(1,in.size(2)) 
              ones(1,in.size(1))  in.size(2)*ones(1,in.size(1)) 1:in.size(2)       1:in.size(2)]';
    
    allpoints = [points_in ; extrapoints];
          
else
    
        
    % create topology
    neighs=[0 0; 1 0; 0 1 ;1 1];
    topology =[];
    for i=1:2
        
        % vertex 1
        [ix,iy] = ndgrid(1+neighs(i,1):in.size(1)-1+neighs(i,1),1+neighs(i,2):in.size(2)-1+neighs(i,2));
        indices1D1=sub2ind(in.size',ix(:),iy(:))';
        
        % vertex 2
        [ix,iy] = ndgrid(1+neighs(i+1,1):in.size(1)-1+neighs(i+1,1),1+neighs(i+1,2):in.size(2)-1+neighs(i+1,2));
        indices1D2=sub2ind(in.size',ix(:),iy(:))';
        
        % vertex 3
        [ix,iy] = ndgrid(1+neighs(i+2,1):in.size(1)-1+neighs(i+2,1),1+neighs(i+2,2):in.size(2)-1+neighs(i+2,2));
        indices1D3=sub2ind(in.size',ix(:),iy(:))';
        
        if (i==1)
            topology = [topology ; indices1D1(:) indices1D2(:) indices1D3(:) ];
        else
            topology = [topology ; indices1D1(:) indices1D3(:) indices1D2(:) ];
        end
        
        %     [ix, iy] [ix+1, iy] [ix, iy+1]
        %     [ix+1, iy] [ix+1, iy+1] [ix, iy+1]        
    end
    
    borders = [1:in.size(1)       1:in.size(1) ones(1,in.size(2))  in.size(1)*ones(1,in.size(2)) 
              ones(1,in.size(1))  in.size(2)*ones(1,in.size(1)) 1:in.size(2)       1:in.size(2)]';
    
    
end

outm = MeshType();
outm.points = allpoints;
outm.npoints = size(outm.points,1);
outm.points = [outm.points zeros(outm.npoints,1)];
outm.triangles = topology;
outm.ntriangles =size(outm.triangles,1 );

% copy the values from the image (using linear interp)

values = in.GetValue(outm.points(:,1:2)',interpolation);

at = AttributeType();
at.attribute='scalars';
at.name='intensity';
at.numComp=1; % between 1 and 4
at.type='float';
at.nelements=outm.npoints; %number of tuples
at.attribute_array = values;
outm.attributes(1)=at;

isBorder=zeros(outm.npoints,1);
indices1Dborder=sub2ind(in.size',borders(:,1),borders(:,2))';
isBorder(indices1Dborder)=1;

at = AttributeType();
at.attribute='scalars';
at.name='border';
at.numComp=1; % between 1 and 4
at.type='float';
at.nelements=outm.npoints; %number of tuples
at.attribute_array = isBorder;
outm.attributes(2)=at;

end

