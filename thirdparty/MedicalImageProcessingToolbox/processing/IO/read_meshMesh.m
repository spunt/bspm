function mesh =read_meshMesh(filename)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType

%from vtkCellType.h
keywordList={'MeshVersionFormatted','Dimension','Vertices','Triangles','Quadrilaterals','Tetrahedra',...
    'Hexahedra','Normals','Tangents','Corners','Edges','Ridges','End'};


if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.vtk', 'Read vtk-file');
    filename = [pathname filename];
end

fid=fopen(filename,'rb');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end

mesh=MeshType();
isEnd=false;

dimension=-1;
nvertices=-1;
vertices=[];
ntriangles=-1;
triangles=[];


while ~isEnd
    kw = nextKeyword(fid,keywordList);
    switch kw{1}
        case keywordList{1}
            % This is currently ignored. states wether is float or double
        case keywordList{2}
            kw = nextKeyword(fid,keywordList);
            dimension = str2num(kw{1});
        case keywordList{3}
            kw = nextKeyword(fid,keywordList);
            nvertices= str2num(kw{1});
            vertices=fscanf(fid,'%f %f %f %f',[dimension+1 nvertices ]);
            vertices = vertices(1:3,:)';
            mesh.npoints = nvertices;
            mesh.points=vertices;
        case keywordList{4}
            kw = nextKeyword(fid,keywordList);
            ntriangles= str2num(kw{1});
            triangles=fscanf(fid,'%f %f %f %f',[dimension+1 ntriangles]);
            triangles = triangles(1:3,:)';
            mesh.ntriangles = ntriangles;
            mesh.triangles=triangles;
        case keywordList{5}
            % quadrilaterals
            kw = nextKeyword(fid,keywordList);
            ntriangles= str2num(kw{1});
            quadrilaterals=fscanf(fid,'%d %d %d %d %d',[dimension+2 ntriangles]);
            quadrilaterals= quadrilaterals(1:4,:)';
            % convert each cuadrilateral into 2 triangles
            triangles = [quadrilaterals(:,1:3) ; quadrilaterals(:,[3 4 1])];
            mesh.ntriangles = size(triangles,1);
            mesh.triangles=triangles;
            
        case keywordList{end}
            %disp(['Finished  '])
            isEnd=true;
        otherwise
            %disp(['Unknown field ' kw{1}])
            
    end
    
end
end


function kw = nextKeyword(fid,keywordList)
str = fgetl(fid);
if strcmp(str,keywordList{end})
    kw=keywordList{end};
    return;
end
[a1,a2]=strtok(str);
while strcmp(a1,'#')
    str = fgetl(fid);
    [a1,a2]=strtok(str);
end

kw{1}=a1;
kw{2}=a2;


end
