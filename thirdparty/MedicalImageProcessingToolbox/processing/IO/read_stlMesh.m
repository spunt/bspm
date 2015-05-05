function mesh =read_stlMesh(filename)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType

if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.stl', 'Read stl-file');
    filename = [pathname filename];
end

fid = fopen(filename,'r','l');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end

% read header
header = fread(fid,80,'uint8');
n_triangular_facets = fread(fid,1,'uint32');

triangles = zeros(n_triangular_facets,3);
normals = zeros(n_triangular_facets,3);
points = [];

for i=1:n_triangular_facets 
    normal = fread(fid,3,'float32');
    v1 = fread(fid,3,'float32');
    v2 = fread(fid,3,'float32');
    v3 = fread(fid,3,'float32');
    abc = fread(fid,1,'uint16'); % attribute byte count
    
    
    % check if the point already exists
    triangles(i,:) = [size(points,1)+1 size(points,1)+2 size(points,1)+3];
    i1=[];
    i2=[];
    i3=[];
    if size(points,1)>0
        dist = (points - ones(size(points,1),1)*v1(:)').^2;
        dist = sum(dist,2);
        i1 = find(dist==0);
        if numel(i1)
            triangles(i,1)=i1;
            triangles(i,2:3)=triangles(i,2:3)-1;
        end
        dist = (points - ones(size(points,1),1)*v2(:)').^2;
        dist = sum(dist,2);
        i2 = find(dist==0);
        if numel(i2)
            triangles(i,2)=i2;
            triangles(i,3)=triangles(i,3)-1;
        end
        dist = (points - ones(size(points,1),1)*v3(:)').^2;
        dist = sum(dist,2);
        i3 = find(dist==0);
        if numel(i3)
            triangles(i,3)=i3;
        end
        %triangles(i,find([i1   i2 i3])) =[numel(i1)*i1   numel(i2)*i2 numel(i3)*i3];
        
    end
    
    points = [points
                v1(1)*(find(nnz(mod(numel(i1)+1,2)))) v1(2)*(find(nnz(mod(numel(i1)+1,2)))) v1(3)*(find(nnz(mod(numel(i1)+1,2))))
                v2(1)*(find(nnz(mod(numel(i2)+1,2)))) v2(2)*(find(nnz(mod(numel(i2)+1,2)))) v2(3)*(find(nnz(mod(numel(i2)+1,2))))
                v3(1)*(find(nnz(mod(numel(i3)+1,2)))) v3(2)*(find(nnz(mod(numel(i3)+1,2)))) v3(3)*(find(nnz(mod(numel(i3)+1,2))))];

    normals(i,:)=normal(:)';
end
fclose(fid);

mesh  = MeshType();
mesh.triangles = triangles;
mesh.points = points;
mesh.ntriangles = size(triangles,1);
mesh.npoints= size(points,1);

mesh.removeDuplicatedPoints();


 
    at = AttributeType();
    at.attribute='normals';
    at.name='normalVectorsFaces';
    at.numComp=3; % between 1 and 4
    at.type='float';
    at.nelements=mesh.ntriangles; %number of tuples
    at.attribute_array = normals;
    
    ind = mesh.find_attribute(at.name);
    if (ind > 0)
        mesh.attributes(ind)=at;
    else
        n_attributes = numel(mesh.attributes);
        if n_attributes==1 && ~numel(mesh.attributes(1).name)
            mesh.attributes(1)=at;
        else
            mesh.attributes(n_attributes+1)=at;
        end
            
    end


end

