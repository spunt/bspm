function c = combineMeshes(m1,m2)
% c = combineMeshes(m1,m2)
%
% This function combines 2 intersecting meshes, i.e. computes the
% intersection line and adds the necessary conectivities
%
% m1, m2 and c are of the class MeshType

% The alorithm works as follows:
% for each face in mesh m1, we see if the defined pane intersects any of
% the triangles of the destination mesh, and if so, if the intersection
% points are within the triangle

c=MeshType(m1);
c.triangles = m1.triangles;
c.points = m1.points;

for i=1:m1.ntriangles
    
    points1 = m1.points(m1.triangles(i,:),:); % each point is one row
    
    for j=1:m2.ntriangles
    
        points2 = m2.points(m2.triangles(j,:),:); % each point is one row
        
        [pts triangle_index ]= intersectionTriangleTriangle(points1,points2);
        if (numel(pts)==0)
            continue;
        end
        % In addition add that the two points define a new connectivity!
        % Think that depending on which triangle they come from, it will be
        % different the new connectivity!
        c.points = [ c.points ;pts'];
       
        % Add to the labels field a new values !
        
        newTriangle = [c.npoints+1 m1.triangles(i,1) m1.triangles(i,2); c.npoints+1 m1.triangles(i,2) m1.triangles(i,3); c.npoints+1 m1.triangles(i,3) m1.triangles(i,1)];
        c.triangles = [c.triangles ; newTriangle];
        
        c.npoints = size(c.points,1);
        
    end

end

% plot results




end