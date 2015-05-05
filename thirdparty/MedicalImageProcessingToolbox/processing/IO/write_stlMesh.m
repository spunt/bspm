function write_stlMesh(filename, mesh)
% Function for writing a mesh in a Visualization Toolkit (VTK) format
% 
% mesh  = write_stlMesh(filename);
%
% examples:
%   mesh=write_stlMesh('mesh.stl');
%   mesh would be of the class MeshType


fid = fopen(filename,'w','l');
if(fid<0)
    fprintf('could not open file %s for writing\n',filename);
    return
end

header = zeros(80,1);
fwrite(fid,header,'uint8'); %header
fwrite(fid,mesh.ntriangles,'uint32');

ind = mesh.find_attribute('normalVectorsFaces');
if (ind <= 0)
    mesh.ComputeNormalsToFaces();
    ind = mesh.find_attribute('normalVectorsFaces');
end

for i=1:mesh.ntriangles
    normal = mesh.attributes(ind).attribute_array(i,:);
    fwrite(fid,normal(:),'float32');
    points = mesh.points(mesh.triangles(i,:),:);
    fwrite(fid,points(1,:)','float32');
    fwrite(fid,points(2,:)','float32');
    fwrite(fid,points(3,:)','float32');
    fwrite(fid,0,'uint16'); % attribute byte count
end
fclose(fid);


end
  
