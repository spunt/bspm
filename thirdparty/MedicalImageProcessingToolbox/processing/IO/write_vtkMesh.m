function write_vtkMesh(filename, mesh)
% Function for writing a mesh in a Visualization Toolkit (VTK) format
% 
% mesh  = write_vtkMesh(filename);
%
% examples:
%   mesh=write_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType


%from vtkCellType.h
cell_types_dict={'VTK_VERTEX','VTK_POLY_VERTEX','VTK_LINE','VTK_POLY_LINE','VTK_TRIANGLE','VTK_TRIANGLE_STRIP',...
    'VTK_POLYGON','VTK_PIXEL','VTK_QUAD','VTK_TETRA','VTK_VOXEL','VTK_HEXAHEDRON','VTK_WEDGE','VTK_PYRAMID',...
    'VTK_PARAMETRIC_CURVE','VTK_PARAMETRIC_SURFACE'} ; % 0 = empty cell 1 to 16
% VTK_QUADRATIC_TETRA                  = 24,

fid=fopen(filename,'wb');
if(fid<0)
    fprintf('could not open file %s for writing\n',filename);
    return
end
fclose(fid);


% write header
dlmwrite(filename, [ '# vtk DataFile Version 3.0' ],'delimiter','');
%if size(mesh.triangles,2)==4 % tetrahedra
%    dlmwrite(filename, [ 'vtk output' ],'delimiter','','-append');
%elseif size(mesh.triangles,2)>3

if (isa(mesh,'QuadTetMeshType'))
    dlmwrite(filename, [ 'vtk output' ],'delimiter','','-append');
elseif size(mesh.triangles,2)>3
    dlmwrite(filename, [ 'vtk output' ],'delimiter','','-append');
else
    dlmwrite(filename, [ 'vtk output' ],'delimiter','','-append');
    %dlmwrite(filename, [ 'DataEntityType: Surface mesh' ],'delimiter','','-append');
end

dlmwrite(filename, [ 'ASCII' ],'delimiter','','-append');
% if size(mesh.triangles,2)==4
%     dlmwrite(filename, [ 'DATASET POLYDATA' ],'delimiter','','-append');
% elseif size(mesh.triangles,2)>3
if (isa(mesh,'QuadTetMeshType'))
    dlmwrite(filename, [ 'DATASET UNSTRUCTURED_GRID' ],'delimiter','','-append');
elseif size(mesh.triangles,2)>3
    dlmwrite(filename, [ 'DATASET UNSTRUCTURED_GRID' ],'delimiter','','-append');
else
    %dlmwrite(filename, [ 'DATASET POLYDATA' ],'delimiter','','-append');
     dlmwrite(filename, [ 'DATASET UNSTRUCTURED_GRID' ],'delimiter','','-append');
end


% write points
dlmwrite(filename, [ 'POINTS ' num2str(mesh.npoints) ' float' ],'delimiter','','-append');
dlmwrite(filename, mesh.points,'delimiter',' ','-append','precision','%.9f');

if (isa(mesh,'QuadTetMeshType'))
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELL_TYPES ' num2str(mesh.ncells) ],'delimiter','','-append');
    dlmwrite(filename, [ ones(mesh.ncells,1)*24 ] ,'delimiter',' ','-append');
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELLS ' num2str(mesh.ncells) ' ' num2str((size(mesh.cells,2)+1)*mesh.ncells)  ],'delimiter','','-append');
    dlmwrite(filename, [ ones(mesh.ncells,1)*size(mesh.cells,2) mesh.cells-1] ,'delimiter',' ','-append','precision', 10);
    
elseif size(mesh.triangles,2)==4 %write tetrahedra
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELL_TYPES ' num2str(mesh.ntriangles) ],'delimiter','','-append');
     if size(mesh.triangles,2)==4
        dlmwrite(filename, [ ones(mesh.ntriangles,1)*10 ] ,'delimiter',' ','-append');
    end
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELLS ' num2str(mesh.ntriangles) ' ' num2str((size(mesh.triangles,2)+1)*mesh.ntriangles)  ],'delimiter','','-append');
    dlmwrite(filename, [ ones(mesh.ntriangles,1)*size(mesh.triangles,2) mesh.triangles-1] ,'delimiter',' ','-append','precision', 10);
%write triangles
elseif size(mesh.triangles,2)>3
    
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELL_TYPES ' num2str(mesh.ntriangles) ],'delimiter','','-append');
    if size(mesh.triangles,2)==4
        dlmwrite(filename, [ ones(mesh.ntriangles,1)*10 ] ,'delimiter',' ','-append');
    end
    
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELLS ' num2str(mesh.ntriangles) ' ' num2str((size(mesh.triangles,2)+1)*mesh.ntriangles)  ],'delimiter','','-append');
    dlmwrite(filename, [ ones(mesh.ntriangles,1)*size(mesh.triangles,2) mesh.triangles-1] ,'delimiter',' ','-append','precision', 10);
else
%     dlmwrite(filename, ' ','delimiter','\n','-append');
%     dlmwrite(filename, [ 'POLYGONS ' num2str(mesh.ntriangles) ' ' num2str(4*mesh.ntriangles)  ],'delimiter','','-append');
%     dlmwrite(filename, [ ones(mesh.ntriangles,1)*3 mesh.triangles-1] ,'delimiter',' ','-append');
%     dlmwrite(filename, ' ','delimiter','\n','-append');
%     dlmwrite(filename, [ 'CELL_DATA ' num2str(mesh.ntriangles)   ],'delimiter','','-append');
%     dlmwrite(filename, [ 'SCALARS scalars short'],'delimiter','','-append');
%     dlmwrite(filename, [ 'LOOKUP_TABLE default '   ],'delimiter','','-append');
%     dlmwrite(filename,  ones(mesh.ntriangles,1) ,'delimiter',' ','-append');
    
    
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELL_TYPES ' num2str(mesh.ntriangles) ],'delimiter','','-append');
    dlmwrite(filename, [ ones(mesh.ntriangles,1)*5 ] ,'delimiter',' ','-append'); % triangles
    dlmwrite(filename, ' ','delimiter','\n','-append');
    dlmwrite(filename, [ 'CELLS ' num2str(mesh.ntriangles) ' ' num2str((size(mesh.triangles,2)+1)*mesh.ntriangles)  ],'delimiter','','-append');
    dlmwrite(filename, [ ones(mesh.ntriangles,1)*size(mesh.triangles,2) mesh.triangles-1] ,'delimiter',' ','-append','precision', 10);
    
    
end

% write attributes

nattributes = numel(mesh.attributes);
total_nfielddatas = 0;
nfielddatas = 0;
if (nattributes==1 && (numel(mesh.attributes.nelements) && mesh.attributes.nelements>0) || nattributes>1)
    dlmwrite(filename, [ ' ' ],'delimiter','','-append');
    dlmwrite(filename, [ 'POINT_DATA ' num2str(mesh.attributes(1).nelements) ],'delimiter','','-append');
    for j=1:nattributes
        if ~numel(lower(mesh.attributes(j).attribute))
            continue;
        end
        switch (lower(mesh.attributes(j).attribute))
            case 'scalars'
                dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' ' mesh.attributes(j).name ' ' mesh.attributes(j).type ],'delimiter','','-append');
                dlmwrite(filename, [ 'LOOKUP_TABLE ' mesh.attributes(j).lookup_table ],'delimiter','','-append');
                dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter','','-append','precision','%.9f');
            case 'field'
                nfielddatas = nfielddatas+1;
                if nfielddatas==1
                    for jj=j:nattributes
                        if (strcmp(lower(mesh.attributes(jj).attribute),'field'))
                            total_nfielddatas=total_nfielddatas+1;
                        end
                    end
                    %dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(nfielddatas)  ],'delimiter','','-append');
                    dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(total_nfielddatas)  ],'delimiter','','-append');
                end
                dlmwrite(filename, [ mesh.attributes(j).name ' ' num2str(mesh.attributes(j).numComp) ' ' num2str(mesh.attributes(j).nelements) ' ' mesh.attributes(j).type ],'delimiter','','-append');
                if strcmp(mesh.attributes(j).type,'short')
                    dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append','precision','%d');    
                else
                    dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append','precision','%.9f');
                end
             case 'vectors'
                 %nfielddatas = nfielddatas+1;
                 dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(nfielddatas) ],'delimiter','','-append');
                 %dlmwrite(filename, [ mesh.attributes(j).name   ' '  num2str(mesh.attributes(j).numComp) ' '  num2str(mesh.attributes(j).nelements) ' '      mesh.attributes(j).type ],'delimiter','','-append');
                 %dlmwrite(filename, [ upper(mesh.attributes(j).attribute) '  ' mesh.attributes(j).name   ' ' mesh.attributes(j).type ],'delimiter','','-append');
                 dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append','precision','%.9f');
            otherwise
                % do nothing
        end


    end

end
  
