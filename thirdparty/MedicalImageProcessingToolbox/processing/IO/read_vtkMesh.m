function mesh =read_vtkMesh(filename)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType

%from vtkCellType.h
cell_types_dict={'VTK_VERTEX','VTK_POLY_VERTEX','VTK_LINE','VTK_POLY_LINE','VTK_TRIANGLE','VTK_TRIANGLE_STRIP',...
    'VTK_POLYGON','VTK_PIXEL','VTK_QUAD','VTK_TETRA','VTK_VOXEL','VTK_HEXAHEDRON','VTK_WEDGE','VTK_PYRAMID',...
    'VTK_PARAMETRIC_CURVE','VTK_PARAMETRIC_SURFACE'} ; % 0 = empty cell


if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.vtk', 'Read vtk-file');
    filename = [pathname filename];
end

fid=fopen(filename,'rb');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end



str = fgetl(fid);
info.Filename=filename;
info.Format=str(3:5); % Must be VTK
info.Version=str(end-2:end);
info.Header = fgetl(fid);
info.DatasetFormat= lower(fgetl(fid));
str = lower(fgetl(fid));
info.DatasetType = str(9:end);

if ~strcmp(lower(info.DatasetType),'polydata')
    
    if strcmp(lower(info.DatasetType),'unstructured_grid')
        disp('Unstructured grid')
    else
        disp('ERROR: Unrecognized vtk data type');
        return
    end
end

mesh = MeshType(0,0);

n_attributes = 0;

while ~feof(fid)
    str=lower(fgetl(fid));
    % 1. Read Geometry/topology ***********************
    % Point data --------------------------------
    if (strfind(str,'points'))
        [tmp1, tmp2] = strtok(str);
        [tmp1, tmp2] = strtok(tmp2);
        npoints = str2num(tmp1);
        mesh.npoints = npoints;
        %disp(['Read ' num2str(npoints) ' points'])
        p = zeros(npoints,3);
        
        [num]=fscanf(fid,'%f');
        index=1;
        for i = 1:npoints
            p(i,1)=   num(index);
            index=index+1;
            p(i,2)=   num(index);
            index=index+1;
            p(i,3)=   num(index);
            index=index+1;
        end
        mesh.points = p;
        clear p;
        
        % Triangle data --------------------------------
    elseif (strfind(str,'polygons'))
        [~, tmp2] = strtok(str);
        [~, nelems] = strtok(tmp2);
        %npolygons = str2double(npolygons);
        nelems = str2double(nelems);
        
        %disp(['Read ' num2str(ntriangles) ' triangles'])
        %mesh.ntriangles = ntriangles;
        %t = zeros(ntriangles,3);
        
        [num, ~]=fscanf(fid,'%d');
        remaining = nelems;
        index=1;
        triangles = [];
        
        while remaining>0
            t=[];
            n = num(index);
            if n==3
                t = num(index+1:index+3)';
                index = index+n+1;
                remaining = remaining-(n+1);
            elseif n==4
                
                t = [num([index+1 index+2 index+3])' ; num([index+3 index+4 index+1])'];
                index = index+n+1;
                remaining = remaining-(n+1);
            else
                disp([ 'Unknown type of cell with ' num2str(n) ' points, remaining ' num2str(remaining)])
                index = index+n+1;
                remaining = remaining-(n+1);
            end
            triangles = [triangles ; t];
            
        end
        
        
        
        
        mesh.triangles = triangles+1; % in matlab, first index is 1;
        mesh.ntriangles = size(mesh.triangles,1);
        clear t;
        % Cell data instead of triangles data --------------------------------
    elseif (strfind(str,'cells'))
        [tmp1, tmp2] = strtok(str);
        [tmp1, tmp2] = strtok(tmp2);
        
        ncells = str2num(tmp1);
        
        %disp(['Read ' num2str(ntriangles) ' triangles'])
        mesh.ntriangles = ncells;
        
        
        [num count]=fscanf(fid,'%d');
        t = zeros(ncells, num(1));
        index=1;
        for i = 1:ncells
            nnodes = num(index);
            if nnodes ~= num(1)
                %disp([ 'Found cell with ' num2str(nnodes ) ' elements, which I am ignoring' ])
                t(i,:)=   -1*ones(1,num(1));
                continue;
            end
            nn=1;
            for nn=1:nnodes
                index=index+1;
                t(i,nn)=   num(index);
            end
            index=index+1;
        end
        
        t(find(t(:,1)==-1),:)=[]; % remove cells which have a different number of nodes
        mesh.triangles = t+1; % in matlab, first index is 1;
        mesh.ntriangles = size(t,1);
        clear t;
        % 2. Read Dataset attributes ********************
        % Point data --------------------------------//////////////
    elseif (strfind(str,'point_data'))
        [tmp1, tmp2] = strtok(str);
        [tmp1, tmp2] = strtok(tmp2);
        npoint_data = str2num(tmp1);
        %disp(['Read  ' num2str(npoint_data) ' scalar point data'])
        
        still_in_point_data = true;
        while(still_in_point_data)
            % build attribute header ---------------
            % read first line
            str=lower(fgetl(fid));   [attribute, tmp2] = strtok(str);
            
            if (strcmp(attribute,'scalars'))
                n_attributes = n_attributes+1;
                at = AttributeType();
                
                [attribute_name, tmp2] = strtok(tmp2);      [type, nc] = strtok(tmp2);
                
                at.attribute = attribute;
                at.name = attribute_name;
                at.type = type;
                if (numel(nc)==0)
                    at.numComp = 1; % default value
                else
                    at.numComp = str2num(nc); % default value
                end
                at.nelements = npoint_data;
                %read second line
                str=lower(fgetl(fid));% LOOKUP_TABLE default
                [lt, name] = strtok(str);
                at.lookup_table=name;
                % fill attribute body ---------------------
                if strcmp(at.type,'short')
                    [num count]=fscanf(fid,'%d');
                elseif (strcmp(at.type,'float') || strcmp(at.type,'double'))
                    [num count]=fscanf(fid,'%f');
                end
                at.attribute_array = zeros(at.nelements,at.numComp);
                index=1;
                for ii = 1:at.nelements
                    for jj=1:at.numComp
                        at.attribute_array(ii,jj)=   num(index);
                        index = index+1;
                    end
                end
                attributes(n_attributes)=at;
                
                
            elseif (strcmp(attribute,'color_scalars'))
                % do nothing
            elseif (strcmp(attribute,'lookup_table'))
                % do nothing
            elseif (strcmp(attribute,'vectors'))
               % continue;
                [tmp1, tmp2] = strtok(str);
                [attribute_name, tmp2] = strtok(tmp2);
                
                n_attributes = n_attributes+1;
                at = AttributeType();
                [type, nc] = strtok(tmp2);
                at.nelements = mesh.npoints;
                
                
                v = zeros(mesh.npoints,3);
                
                [num]=fscanf(fid,'%f');
                index=1;
                for i = 1:mesh.npoints
                    v(i,1)=   num(index);
                    index=index+1;
                    v(i,2)=   num(index);
                    index=index+1;
                    v(i,3)=   num(index);
                    index=index+1;
                end
                
                at.attribute_array = v;
                at.type = type;
                at.numComp = 3;
                at.name = attribute_name;
                at.attribute = attribute;
                attributes(n_attributes)=at;
                clear v;
            elseif (strcmp(attribute,'normals'))
                
            elseif (strcmp(attribute,'tensors'))
                % do nothing
            elseif (strcmp(attribute,'texture_coordinates'))
                % do nothing
            elseif (strcmp(attribute,'field'))
                
                [fieldname, numArrays] = strtok(tmp2);
                
                for i=1:str2num(numArrays)
                    n_attributes = n_attributes+1;
                    at = AttributeType();
                    at.attribute = attribute;
                    
                    %read second line
                    str=lower(fgetl(fid));
                     if str<0
                         continue;
                     end
                    C = textscan(str, '%s %d %d %s');
                    at.name =C{1}{1};
                    at.numComp=C{2};
                    at.nelements=C{3};
                    at.type =C{4}{1};
                    
                    
                    % fill attribute body ---------------------
                    if strcmp(at.type,'short')
                        [num count]=fscanf(fid,'%d');
                    elseif (strcmp(at.type,'float') || strcmp(at.type,'double'))
                        [num count]=fscanf(fid,'%f');
                    end
                    at.attribute_array = zeros(at.nelements,at.numComp);
                    index=1;
                    for ii = 1:at.nelements
                        for jj=1:at.numComp
                            at.attribute_array(ii,jj)=   num(index);
                            index = index+1;
                        end
                    end
                    attributes(n_attributes)=at;
                    %disp(['Read  ' num2str(at.nelements) ' field data'])
                end
            else
                still_in_point_data = false;
            end
        end
        
        % end of reading point data ----------------------/////////
    elseif (strfind(str,'triangle_strips'))
        
        [~, tmp2] = strtok(str);
        [~, nelems] = strtok(tmp2);
        nelems = str2double(nelems);
        
        [num, ~]=fscanf(fid,'%d');
        remaining = nelems;
        index=1;
        triangles = [];
        while remaining>0
            n = num(index);
            
            t = zeros(n-2,3);
            for k=1:2:(n-2)
                t(k,:)=num([index+k index+k+1 index+k+2])';
            end
            for k=2:2:(n-2)
                t(k,:)=num([index+k index+k+2 index+k+1])';
            end
            index = index+n+1;
            remaining = remaining-(n+1);
            triangles = [triangles ; t];
        end
        
        mesh.triangles = triangles+1; % in matlab, first index is 1;
        mesh.ntriangles = size(mesh.triangles,1);
        clear t;
        
    end
    
    
end
if (n_attributes>0)
    mesh.attributes = attributes;
end
fclose(fid);
end

