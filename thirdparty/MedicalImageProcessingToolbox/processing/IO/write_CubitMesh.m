function write_CubitMesh(filenameFac,m)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType





dlmwrite(filenameFac,[m.npoints m.ntriangles],'delimiter','\t')
dlmwrite(filenameFac,[1:m.npoints ;  m.points']','-append','delimiter','\t');

%dlmwrite(filenameFac,m.ntriangles,'-append','delimiter','\t')
dlmwrite(filenameFac,[1:m.ntriangles ;  m.triangles']','-append','delimiter','\t');

end

