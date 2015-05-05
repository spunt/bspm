function d = distanceMeshToMesh( m1,m2 )
%DISTANCEMESHTOMESH calculate distance between meshes
% works best if m2 < m1 in size

% case one: mean shortest distance: for each point of m2, we pick the
% closest point in m1

distances = zeros(m2.npoints,1);

for i=1:m2.npoints
    
   dist = m1.points - ones(m1.npoints,1)*m2.points(i,:);
   norm2 = sqrt(dot(dist,dist,2));
   distances(i) = min(sqrt(norm2));
    
end

d = mean(distances);


end

