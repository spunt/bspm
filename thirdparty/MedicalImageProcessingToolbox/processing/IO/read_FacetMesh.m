function m =read_FacetMesh(filename)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType


fid = fopen(filename,'r','l');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end

m = MeshType();

%read coordinates
npoints = str2num(fgetl(fid));
tline = fscanf(fid,'%d\t%f\t%f\t%f\n', 4*npoints);
m.points= reshape(tline,4,npoints)';
m.points = m.points(:,2:4);

%read facets
ntriangles = str2num(fgetl(fid));
tline = fscanf(fid,'%d\t%d\t%d\t%d\n', 4*ntriangles);
m.triangles= reshape(tline,4,ntriangles)';
m.triangles = m.triangles(:,2:4);

fclose(fid);
m.npoints = size(m.points,1);
m.ntriangles= size(m.triangles,1);


end

