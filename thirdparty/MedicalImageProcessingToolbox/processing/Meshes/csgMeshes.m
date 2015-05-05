function output = csgMeshes(m1,m2,varargin)
% outMesh = csgMeshes(m1,m2)
%   
%   Computes different operations of constructive solid geometry between
%   two meshes: clip, Union, difference and intersection
%
% m1, m2 and output are of the class MeshType

operation = 'intersection';
dbg = false;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'op'))
        operation= varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'debug'))
        dbg = true;
    end
end

c=MeshType(m1);

% Compute, for each facet of the input mesh with less facets, the 
% AABB-tree (Axis aligned Bounding Box tree)

if (m2.ntriangles < m1.ntriangles)
    meshSmall = m2;
    meshBig = m1;
else
    meshSmall = m1;
    meshBig = m2;
end

treeSmall = AABB_tree(meshSmall);
treeSmall.BuildTree();

treeBig = AABB_tree(meshBig);
treeBig.BuildTree();

% Using the trees, compute the intersections

% start up with the larger scale and move deeper into the tree. if the 
% AABBs fail to overlap in any direction, then they must not intersect at 
% all; if the AABBs overlap in all directions, then they must intersect.


treeSmall.intersectTrees(treeBig);
% pointsx = treeSmall.intersection_path(:,1);
% pointsy = treeSmall.intersection_path(:,2);
% pointsz = treeSmall.intersection_path(:,3);

% viewMesh(treeSmall.rootNode.mesh,'wireframeSurface');
% viewMesh(treeBig.rootNode.mesh,'color',[0 0 1],'wireframeSurface');
% hold on;
% line(pointsx(treeSmall.newConectivity'),...
%     pointsy(treeSmall.newConectivity'),...
%     pointsz(treeSmall.newConectivity'),'LineWidth',4,'color',[0 1 0]);
% hold off;

outMesh = treeSmall.generateOutputMesh(operation);

%figure; viewMesh(outMesh{1},'wireframeSurface','color',[0 0 1])
%hold on;  viewMesh(outMesh{2},'wireframeSurface','color',[0 1 0]); hold off;
%hold on;  viewMesh(outMesh{3},'wireframeSurface','color',[1 0 0]); hold off;
%hold on;  viewMesh(outMesh{4},'wireframeSurface','color',[0.5 0.5 0]); hold off;
%alpha(0.7);

output = MeshType(outMesh{1});
if strcmp(operation,'intersection') || strcmp(operation,'difference') || strcmp(operation,'union')
    output.triangles = [output.triangles ; outMesh{3}.triangles+output.npoints];
    output.points = [output.points ; outMesh{3}.points];
    output.ntriangles = size(output.triangles,1);
    output.npoints= size(output.points,1);
    output.removeDuplicatedPoints();
end
% TODO: build the appropriate output:
% If union, difference or intersection -> outMesh{1} + outMesh{3}
% if clip, -> outMesh 1



end