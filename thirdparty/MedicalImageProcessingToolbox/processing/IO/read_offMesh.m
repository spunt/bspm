function mesh =read_offMesh(filename)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType

if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.off', 'Read off-file');
    filename = [pathname filename];
end

fid=fopen(filename,'rb');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end



str = fgetl(fid);
tmp=fscanf(fid,'%d %d %d',3);
npoints = tmp(1);
ntriangles= tmp(2);
mesh = MeshType(npoints,ntriangles);
points=fscanf(fid,'%f',3*npoints);
mesh.points = reshape(points,3,[])';
triangles=fscanf(fid,'%f',4*ntriangles);
triangles = reshape(triangles,4,[])';
mesh.triangles = triangles(:,2:end)+1;

fclose(fid);
end

