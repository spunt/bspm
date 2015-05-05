function out = cylinderMesh(axis,centre,radius, height,resolution)

angles = (0:(resolution-1))/(resolution-1)*2*pi;

points = [radius*cos(angles) ; radius*sin(angles); zeros(size(angles))];
pointsUp = points + [0 0 height/2]'*ones(1,numel(angles));
pointsDown = points - [0 0 height/2]'*ones(1,numel(angles));

[x,y] = vtkMathPerpendiculars(axis,pi/2);

M = [x y axis centre; 0 0 0 1];

points3D = M*[pointsUp  pointsDown ; ones(1,2*numel(angles))];
    
out = MeshType();
out.npoints = size(points3D,2 );
out.points = points3D(1:3,:)';

% topology
triangles = [];
for i=1:numel(angles)-1
    triangles  = [triangles
        i i+numel(angles) i+1;
        i+1  i+numel(angles) i+numel(angles)+1
        ];
end

out.triangles = triangles;
out.ntriangles = size(triangles,1);

out.attributes.attribute = 'scalars';
out.attributes.name = 'color';
out.attributes.nelements = out.npoints;
out.attributes.attribute_array = zeros(out.npoints,1);
end
