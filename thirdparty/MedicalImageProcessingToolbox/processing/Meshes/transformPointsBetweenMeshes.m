function transformed_points = transformPointsBetweenMeshes( original_mesh, transformed_mesh, points )
%TRANSFORMPOINTSBETWEENMESHES Transforms points on the surface of a 2D mesh
% to points on the surface of another mesh (2D or 3D) that is conformal
% with the first one.
% each row of 'points' is a point
%
% points must be provided as row vectors


bt = original_mesh.ComputeBounds();
extent = bt(2:2:end)-bt(1:2:end);
nodim = find(extent==0);
availabledims = setdiff(1:numel(extent),nodim);

DT_original = delaunayTriangulation(original_mesh.points(:,availabledims));

% for debug
if (0)
    figure,
    triplot(DT_original)
    axis equal
    IC = incenter(DT_original);
    hold on
    numtri = size(DT_original,1);
    trilabels = arrayfun(@(P) {sprintf('T%d', P)}, (1:numtri)');
    Htl = text(IC(:,1),IC(:,2),trilabels,'FontWeight','bold', ...
        'HorizontalAlignment','center','Color','blue');
    plotpoints2(points(:,:),'r.');
    hold off
end

bt = transformed_mesh.ComputeBounds();
extent = bt(2:2:end)-bt(1:2:end);
nodim = find(extent==0);
availabledims2 = setdiff(1:numel(extent),nodim);

%DT_transformed = delaunayTriangulation(transformed_mesh.points(:,availabledims),DT_original.edges);



transformed_points = NaN(size(points,1),numel(availabledims2));

triangle_id = DT_original.pointLocation(points(:,availabledims));
valid_points = find(triangle_id==triangle_id);
barycentric_coords = DT_original.cartesianToBarycentric(triangle_id(valid_points),points(valid_points,availabledims));

vertices = DT_original.ConnectivityList(triangle_id(valid_points),:);
a3D = transformed_mesh.points(vertices(:,1),availabledims2);
b3D = transformed_mesh.points(vertices(:,2),availabledims2);
c3D = transformed_mesh.points(vertices(:,3),availabledims2);
transformed_points(valid_points,:)=a3D.*(barycentric_coords(:,1)*ones(1,numel(availabledims2)))+...
    b3D.*(barycentric_coords(:,2)*ones(1,numel(availabledims2)))+...
    c3D.*(barycentric_coords(:,3)*ones(1,numel(availabledims2)));

end

