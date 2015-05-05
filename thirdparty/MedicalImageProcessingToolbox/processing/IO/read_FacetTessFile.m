function m =read_FacetTessFile(filenameFacets,filenameTESS)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType



points = importdata(filenameFacets,' ',1);
%points.data
triangles= importdata(filenameTESS,' ',1);

m = MeshType();
m.points = points.data;
m.triangles= triangles.data;
m.npoints = size(points.data,1);
m.ntriangles= size(triangles.data,1);

end

