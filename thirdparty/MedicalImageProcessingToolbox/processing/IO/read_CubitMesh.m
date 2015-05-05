function m =read_CubitMesh(filenameFac)
% Function for reading a mesh in a Cubit *.fac format
%
% mesh  = read_vtkMesh(filenameFac);
%
% examples:
%   mesh=read_vtkMesh('volume.fac');
%   mesh would be of the class MeshType

opt.fmt='%f';
opt.header = [1 2];
opt.n=4;

   f = fopen(filenameFac);
   n = fscanf(f,'%i',opt.header);
   data = fscanf(f,opt.fmt,[opt.n Inf]);
   data = data';
   fclose(f);





m = MeshType();
m.points = data(1:n(1),2:4);
m.triangles= data(n(1)+1:end,2:4);
m.npoints = size(m.points,1);
m.ntriangles= size(m.triangles,1);

end

