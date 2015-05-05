function  outmesh = intersectionSpheremeshPlane(m, plane)
%INTERSECTIONSPHEREMESHPLANE Intersect a mesh with a plane.
%
%   SLICE = polyhedronSlice(NODES, FACES, PLANE)
%   M is a mesh of MeshType
%   PLANE: a plane representation [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2].
%   return the intersection polygon of the polyhedra with the plane, in the
%   form of a set of ordered points.
%
%   Works only for convex polyhedra.
%
%   Example
%   polyhedronSlice
%
%   See also
%   polyhedra, clipConvexPolyhedronHP
%
% ------
% Author: Alberto Gomez, inspired by the work by 
%   [ David Legland
%   e-mail: david.legland@nantes.inra.fr
%   Created: 2007-09-18,    using Matlab 7.4.0.287 (R2007a)
%   Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
%   ]


% compute edges of the sphere
inds = zeros(3*m.ntriangles,2); % indices of edges
triangs = reshape(repmat(1:m.ntriangles,3,1),1,[])'; % indices of the corresponding triangle
edgeOrder = repmat(1:3,1,m.ntriangles)'; % indices of the corresponding triangle
for f=1:m.ntriangles
    face =  m.triangles(f,:)';
    current = 3*(f-1)+1;
    %inds(current:(current+2),:) =  sort([face face([2:end 1])], 2);
    inds(current:(current+2),:) =  [face face([2:end 1])];
    clear face;
end

[nonRepeatedEdges i_1 i_2]  = unique(inds, 'rows'); % remove duplicated rows (only for points; triangles need duplicated rows)

edges.p1 = m.points(nonRepeatedEdges(:,1), :)';
edges.p2 = m.points(nonRepeatedEdges(:,2), :)';

% intersection of edges with plane
[points indices] = intersectEdgePlane(edges, plane); % one point per edge
repeatedIndices = indices(i_2);
repeatedPoints = points(:,i_2);

repeatedPoints = repeatedPoints(:,repeatedIndices);
triangs = triangs(repeatedIndices);
edgeOrder = edgeOrder(repeatedIndices);
% Build the new topology


outmesh=MeshType(m);
outmesh.triangles=m.triangles;
outmesh.points=m.points;

edgesMatrix = [1 2; 2 3; 3 1];

index_of_border_vertices=[];

for i = unique(triangs)'
   
    old_triangle = outmesh.triangles(i,:);
    outmesh.triangles(i,:)=[NaN NaN NaN]; % remove previous triangle
    
    % points of the new edge
    
    p_index = find(triangs==i);
    if (numel(p_index) <2)
        continue;
    end
    p1 = repeatedPoints(:,p_index(1));
    p2 = repeatedPoints(:,p_index(2));
    
    
    % project into 2D
    corners = outmesh.points(old_triangle,:);
    normal = cross(corners(2,:)-corners(1,:), corners(3,:)-corners(1,:))';
    normal = normal/norm(normal);
    
    [x_ y_] = vtkMathPerpendiculars(normal, pi/2);
    
    M = [x_ y_ normal];
    
    mypoints2D = M\ [outmesh.points(m.triangles(i,[1 2 3]),:); p1';p2']'; 
    % make sure it is convex
    centre2D = mean(mypoints2D ,2);
    mypoints2D_tmp = mypoints2D -centre2D*ones(1,size(mypoints2D,2));
    mypoints2D_tmp = mypoints2D_tmp*0.6;
    mypoints2D(:,1:3) = mypoints2D_tmp(:,1:3)+centre2D*ones(1,3);
    clear mypoints2D_tmp;
     edgeNotCut = setxor([1 2 3], edgeOrder(p_index) );
    verticesNotCut = old_triangle(edgesMatrix(edgeNotCut,:));
    oppositVertex = setxor([1 2 3],edgesMatrix(edgeNotCut,:));
    
    DT = DelaunayTri(mypoints2D(1:2,:)');
    DT.Constraints = [ 4 5];
    
    
    old_n_points = outmesh.npoints;    
    ps = [old_triangle old_n_points+1  old_n_points+2];
    newTriangles = ps(DT.Triangulation);    
    outmesh.points = [ outmesh.points; p1' ; p2'];               
    outmesh.triangles = [ outmesh.triangles; newTriangles];
    outmesh.npoints = size(outmesh.points,1);
    
    
   % test
   % p=outmesh.points(m.triangles(i,[1 2 3 1]),:); plot3(p(:,1),p(:,2),p(:,3),'k*-');
   % hold on; p=outmesh.points(newTriangles(1,[1 2 3 1]),:); plot3(p(:,1),p(:,2),p(:,3),'ro-','LineWidth',2); hold off;
   % hold on; p=outmesh.points(newTriangles(2,[1 2 3 1]),:); plot3(p(:,1),p(:,2),p(:,3),'go-','LineWidth',2); hold off;
   % hold on; p=outmesh.points(newTriangles(3,[1 2 3 1]),:); plot3(p(:,1),p(:,2),p(:,3),'bo-','LineWidth',2); hold off;
end


index_of_border_vertices = (m.npoints+1):outmesh.npoints;

circleAttribute = zeros(size(outmesh.npoints));
circleAttribute(index_of_border_vertices)=1;

 outmesh.addVertexAttribute( circleAttribute,'intersection');


% Remove old triangles
outmesh.triangles(outmesh.triangles(:,1)~=outmesh.triangles(:,1),:)=[];
outmesh.ntriangles= size(outmesh.triangles,1);

outmesh.npoints = size(outmesh.points,1);









end
