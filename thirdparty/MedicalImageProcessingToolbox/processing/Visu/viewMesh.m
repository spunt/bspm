function viewMesh(mesh,varargin)
% VIEWMESH mesh display using matlab
%
% viewMesh(mesh);
% viewMesh(mesh, option1,val1, ...)
%
% mesh     =   mesh from the class MeshType or subclasses
% options and values are:
%   'showvectors'  --> val,val2; val = scale factor for vectors , val2
%   number of attribute
%   'labelColor'  --> val = if this is active, mesh will be colored by the labels
%   field 'val'.
%   'noplot'    --> if this is active, the image is not plotted (used to accumulate several actors into a same vtkwindow)
%   'color'    --> [r g b] color to display mesh [0,1]
%   'triangles' --> [] array with the triangles indices that are displayed
%
%   TODO continue working on this class

showvectors = false;
scale = 1;
attributeVector=1;
noplot = false;
doWireFrame=false;
opacity = 1;
labelColor=false;
attributeColor=1; % of all the attributes, the one to be used for labelling
color = [0 0 1];
edgeColor=[0 0 0];
linestyle = 'none';
%linestyle = '-'; % uncomment to show edges
i = 1;
faceColor=[1 0 0];
axes_handle = gca;
change_axes=true;
triangles_indices = 1:mesh.ntriangles;
tag = 'Mesh';
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'showvectors'))
            showvectors = true;
            scale = varargin{i+1};
            attributeVector= varargin{i+2};
            i = i+1;
    elseif(strcmp( varargin{i} , 'labelColor'))
            labelColor = true;
            attributeColor= varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'triangles'))
            triangles_indices= varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'axes'))
            change_axes=false;
            axes_handle= varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'wireframe'))
            doWireFrame=true;
            faceColor='none';
            edgeColor = color;
            linestyle = '-';
    elseif(strcmp( varargin{i} , 'wireframeSurface'))
            linestyle = '-';
            edgeColor=[0 0 0];
    elseif(strcmp( varargin{i} , 'color'))
            color= varargin{i+1};
            faceColor = color;
            i = i+1;
    elseif(strcmp( varargin{i} , 'Tag'))
            tag= varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'opacity'))
            opacity= varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'noplot'))
            noplot = true;
    end
    i = i+1;
end

if (doWireFrame)
    faceColor = 'none';
end
   
   
   hold on
    p=0;
    


    if (labelColor)
        colors = mesh.attributes(attributeColor).attribute_array;
         p = trimesh(mesh.triangles(triangles_indices,:), mesh.points(:,1),mesh.points(:,2),mesh.points(:,3), colors,...
             'FaceColor','interp','EdgeColor',[0 0 0],...
             'FaceLighting','gouraud','Parent',axes_handle,'Tag',tag); % 'FaceLighting','phong');
    else
        
         p = trimesh(mesh.triangles(triangles_indices,:), mesh.points(:,1),mesh.points(:,2),mesh.points(:,3),'FaceColor',faceColor,...
             'LineStyle',linestyle,'EdgeColor',edgeColor,'FaceLighting','flat','Parent',axes_handle,'Tag',tag);
        
    end
        
            
        
       
    
    
    if (showvectors)
        % quiver the vectors 
        if size(mesh.attributes(attributeVector).attribute_array,1)==mesh.npoints
            quiver3(mesh.points(:,1), mesh.points(:,2), mesh.points(:,3), mesh.attributes(attributeVector).attribute_array(:,1), mesh.attributes(attributeVector).attribute_array(:,2), mesh.attributes(attributeVector).attribute_array(:,3),scale,'Color',[0 0 0]);
        elseif size(mesh.attributes(attributeVector).attribute_array,1)==mesh.ntriangles
            % find triangle centroids
            pointsOfTriangles = cat(3,mesh.points(mesh.triangles(:,1),:),mesh.points(mesh.triangles(:,2),:),mesh.points(mesh.triangles(:,3),:));
            medicenters = mean(pointsOfTriangles,3);
            quiver3(medicenters(:,1), medicenters(:,2), medicenters(:,3), mesh.attributes(attributeVector).attribute_array(:,1), mesh.attributes(attributeVector).attribute_array(:,2), mesh.attributes(attributeVector).attribute_array(:,3),scale,'Color',[0 0 0]);
        end
    end
    
    
    alpha(p,opacity);
    light('Position',[1 1 0],'Style','infinite'); 
    light('Position',[1 0 1],'Style','infinite'); 
    % axis equal;
    if change_axes
        axis off;
    end
     
     hold off;
 
end