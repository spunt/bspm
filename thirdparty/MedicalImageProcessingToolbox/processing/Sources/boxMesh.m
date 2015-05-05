function out = boxMesh(centre,dimensions)

% dimensions is a 3 element vector with depth, width and height.

[ix,iy,iz]=ndgrid([-.5 .5],[-.5 .5],[-.5 .5]);

points = centre*ones(1,8) +[ix(:) iy(:) iz(:)]'.*(dimensions*ones(1,8));

out = MeshType();
out.npoints = size(points,2 );
out.points = points';

out.triangles = [1 2 3
    3 2 4
    7 6 5
    6 7 8
    1 2 5
    2 6 5
    4 3 7
    7 8 4
    3 1 7
    1 5 7
    2 4 8
    8 6 2];

out.ntriangles = size(out.triangles,1);

out.attributes.attribute = 'scalars';
out.attributes.name = 'color';
out.attributes.nelements = out.npoints;
out.attributes.attribute_array = zeros(out.npoints,1);
end
