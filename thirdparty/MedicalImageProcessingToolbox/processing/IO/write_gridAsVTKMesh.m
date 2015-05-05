function write_gridAsVTKMesh(filename, im)
% Function for writing an image as a grid mesh in a Visualization Toolkit (VTK) format
% 
%  write_gridAsVTKMesh(filename, mesh);
%

%from vtkCellType.h
cell_types_dict={'VTK_VERTEX','VTK_POLY_VERTEX','VTK_LINE','VTK_POLY_LINE','VTK_TRIANGLE','VTK_TRIANGLE_STRIP',...
    'VTK_POLYGON','VTK_PIXEL','VTK_QUAD','VTK_TETRA','VTK_VOXEL','VTK_HEXAHEDRON','VTK_WEDGE','VTK_PYRAMID',...
    'VTK_PARAMETRIC_CURVE','VTK_PARAMETRIC_SURFACE'} ; % 0 = empty cell 1 to 16


fid=fopen(filename,'wb');
if(fid<0)
    fprintf('could not open file %s for writing\n',filename);
    return
end



% write header
dlmwrite(filename, [ '# vtk DataFile Version 3.0' ],'delimiter','');
%if size(mesh.triangles,2)>3
 dlmwrite(filename, [ 'vtk output' ],'delimiter','','-append');
%else
% dlmwrite(filename, [ 'DataEntityType: Surface mesh' ],'delimiter','','-append');
%end
dlmwrite(filename, [ 'ASCII' ],'delimiter','','-append');

%if size(mesh.triangles,2)>3
    dlmwrite(filename,  'DATASET UNSTRUCTURED_GRID' ,'delimiter','','-append');
%else
%    dlmwrite(filename, [ 'DATASET POLYDATA' ],'delimiter','','-append');
%end


% write points
points = im.GetPosition(1:prod(im.size));
dlmwrite(filename, [ 'POINTS ' num2str(prod(im.size)) ' float' ],'delimiter','','-append');
dlmwrite(filename, points','delimiter',' ','-append');


% write lines (calculate topology)

[ixx,iyx,izx]=ndgrid(1:im.size(1)-1,1:im.size(2),1:im.size(3));
[ixy,iyy,izy]=ndgrid(1:im.size(1),1:im.size(2)-1,1:im.size(3));
[ixz,iyz,izz]=ndgrid(1:im.size(1),1:im.size(2),1:im.size(3)-1);

indices0x = [ixx(:) iyx(:) izx(:)];
indices1x = [ixx(:)+1 iyx(:) izx(:)];
indices0y = [ixy(:) iyy(:) izy(:)];
indices1y = [ixy(:) iyy(:)+1 izy(:)];
indices0z = [ixz(:) iyz(:) izz(:)];
indices1z = [ixz(:) iyz(:) izz(:)+1];


indices_0x=sub2ind(im.size',indices0x(:,1),indices0x(:,2),indices0x(:,3));
indices_1x=sub2ind(im.size',indices1x(:,1),indices1x(:,2),indices1x(:,3));
indices_0y=sub2ind(im.size',indices0y(:,1),indices0y(:,2),indices0y(:,3));
indices_1y=sub2ind(im.size',indices1y(:,1),indices1y(:,2),indices1y(:,3));
indices_0z=sub2ind(im.size',indices0z(:,1),indices0z(:,2),indices0z(:,3));
indices_1z=sub2ind(im.size',indices1z(:,1),indices1z(:,2),indices1z(:,3));

topology = [indices_0x(:) indices_1x(:)
            indices_0y(:) indices_1y(:)
            indices_0z(:) indices_1z(:)];


%nlines = prod(im.size-1);
nlines = size(topology,1);
celltype = 3;
% if size(mesh.triangles,2)>3
     
     dlmwrite(filename, ' ','delimiter','\n','-append');
     dlmwrite(filename, [ 'CELL_TYPES ' num2str(nlines) ],'delimiter','','-append');
     
     dlmwrite(filename, [ ones(nlines,1)*celltype ] ,'delimiter',' ','-append');
     
     
     dlmwrite(filename, ' ','delimiter','\n','-append');
     dlmwrite(filename, [ 'CELLS ' num2str(nlines) ' ' num2str((2+1)*nlines)  ],'delimiter','','-append');
     dlmwrite(filename, [ ones(nlines,1)*2 topology-1] ,'delimiter',' ','-append');
% else
%    dlmwrite(filename, [ 'POLYGONS ' num2str(nlines) ' ' num2str(3*nlines)  ],'delimiter','','-append');
%    dlmwrite(filename, [ ones(nlines,1)*2 topology-1] ,'delimiter',' ','-append');
%end

% This I do not need for now

% nattributes = numel(mesh.attributes);
% 
% nfielddatas = 0;
% if (nattributes==1 && mesh.attributes.nelements>0) || nattributes>1
% %if (0)
%     dlmwrite(filename, [ ' ' ],'delimiter','','-append');
%     dlmwrite(filename, [ 'POINT_DATA ' num2str(mesh.attributes(1).nelements) ],'delimiter','','-append');
%     for j=1:nattributes
%         if ~numel(lower(mesh.attributes(j).attribute))
%             continue;
%         end
%         switch (lower(mesh.attributes(j).attribute))
%             case 'scalars'
%                 dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' ' mesh.attributes(j).name ' ' mesh.attributes(j).type ],'delimiter','','-append');
%                 dlmwrite(filename, [ 'LOOKUP_TABLE ' mesh.attributes(j).lookup_table ],'delimiter','','-append');
%                 dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter','','-append');
%             case 'field'
%                 nfielddatas = nfielddatas+1;
%                 dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(nfielddatas)  ],'delimiter','','-append');
%                 dlmwrite(filename, [ mesh.attributes(j).name ' ' num2str(mesh.attributes(j).numComp) ' ' num2str(mesh.attributes(j).nelements) ' ' mesh.attributes(j).type ],'delimiter','','-append');
%                 dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append');
%              case 'vectors'
%                  %nfielddatas = nfielddatas+1;
%                  %dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(nfielddatas) ],'delimiter','','-append');
%                  dlmwrite(filename, [ upper(mesh.attributes(j).attribute) '  ' mesh.attributes(j).name   ' ' mesh.attributes(j).type ],'delimiter','','-append');
%                  dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append');
%             otherwise
%                 % do nothing
%         end
% 
% 
%     end

end
  
