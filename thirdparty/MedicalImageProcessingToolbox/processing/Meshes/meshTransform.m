function meshOut = meshTransform(meshIn, matrix, varargin)
% converts a mesh into a binary image (rasterization), on the grid defined
% 
%   Options:
%   'thickness' , n     ->  n is the number of voxel sizes at each size of
%   the mesh that will be painted (by default, n=1/2)

dbg=false;
invert = false;
i=1;
while (i <= size(varargin,2))
    if(strcmp( varargin{i} , 'debug'))
        dbg= true;
    elseif(strcmp( varargin{i} , 'invert'))
        invert= true;
    end
    i = i+1;
end


meshOut = MeshType(meshIn);
if invert 
    matrix = inv(matrix);
end
meshOut.points = (matrix * [meshIn.points' ; ones(1,meshIn.npoints)])';
meshOut.points = meshOut.points(:,1:3);

end