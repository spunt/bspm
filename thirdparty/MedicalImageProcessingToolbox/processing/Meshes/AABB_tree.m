classdef AABB_tree < handle
    % This class defines a AABB_tree
    % by Alberto Gomez, 2012
    %
    
    properties(GetAccess = 'public', SetAccess = 'public')
        
        rootNode = AABB_node();
        mesh; % structure with input Meshes
        max_tree_depth=Inf; % we want that each ending leaf is a single triangle
        min_vertices=3;
        left_children=0;
        right_children=0;
        curr_max_depth=0;
        dbg=false;
        computedIntersections=[]; % list of all triangle vs triangle intersections that have already been computed
        intersection_path=[];
        newConectivity =[];
        nodes_per_level =zeros(30,1); % very rare that this will be larger than 10 levels!
        
    end
    
    
    
    
    methods(Access = public)
        %constructor
        function obj = AABB_tree(mesh)
            if(nargin > 0)
                obj.mesh{1} = mesh;
                initialDepth=0;
                triangle_identifier = (1:mesh.ntriangles)';
                obj.rootNode = AABB_node(mesh,initialDepth,triangle_identifier,1);
                
            end
        end
        
        function DebugOn(obj)
            obj.dbg = true;
        end
        
        function DebugOff(obj)
            obj.dbg = false;
        end
        
        function BuildTree(obj,varargin)
            
            if (numel(varargin)==0)
                if obj.dbg
                    disp('Processing root')
                end
                currentNode = obj.rootNode;
            else
                if obj.dbg
                    disp([ 'Processing depth ' num2str(obj.curr_max_depth) ])
                end
                currentNode = varargin{1};
            end
            
            if (currentNode.depth+1 > obj.curr_max_depth)
                obj.curr_max_depth=currentNode.depth+1;
            end
            
            if (currentNode.mesh.ntriangles<2 || currentNode.depth>=obj.max_tree_depth)
                return;
            end
            
            
            % Continue computing stuff  -------------------
            
            % Find the longest axii of the node's box
            %int iAxis=Box.GetLongestAxis();
            iAxis = currentNode.FindBestAxis();
            
            
            
            %btw, things that could go wrong- if the mesh doesn't
            %send over the triangles correctly, then you might see
            %huge boxes that are misaligned (bad leaves).
            %things to check is making sure the vertex buffer is
            %correctly aligned along the adjancey buffers, etc
            
            % for all vertices
            
            pointsOfTriangles = cat(3,currentNode.mesh.points(currentNode.mesh.triangles(:,1),:),...
                currentNode.mesh.points(currentNode.mesh.triangles(:,2),:),...
                currentNode.mesh.points(currentNode.mesh.triangles(:,3),:));
            faceCenter = mean(pointsOfTriangles,3); % contains the center of each face
            
            %get center of the box for longest axis
            %float fSplit = ( max[ iAxis ] + min[ iAxis ] ) / 2.f ;
            % center =  (currentNode.bounds([1 3 5])+currentNode.bounds([2 4 6]))/2;
            center =  mean(faceCenter); % center of all the box
            
            
            left = find(faceCenter(:,iAxis) <= center(iAxis));
            right = setdiff((1:size(currentNode.mesh.triangles,1))' , left);
            
            %Build child nodes
            if (numel(left) > 0)
                % all the points in the left side
                childrenMesh = currentNode.mesh.removeTriangle(right);
                obj.left_children = obj.left_children+1;
                triangles_identifier = currentNode.unique_triangles_id(left);
                obj.nodes_per_level(currentNode.depth+1)=obj.nodes_per_level(currentNode.depth+1)+1;
                currentNode.leftNodes = AABB_node(childrenMesh, currentNode.depth+1, triangles_identifier,obj.nodes_per_level(currentNode.depth+1));                
                obj.BuildTree(currentNode.leftNodes);
            end
            if (numel(right) > 0)
                
                childrenMesh = currentNode.mesh.removeTriangle(left);
                obj.right_children =  obj.right_children+1;
                triangles_identifier = currentNode.unique_triangles_id(right);
                obj.nodes_per_level(currentNode.depth+1)=obj.nodes_per_level(currentNode.depth+1)+1;
                currentNode.rightNodes =  AABB_node( childrenMesh, currentNode.depth+1,triangles_identifier,obj.nodes_per_level(currentNode.depth+1));                
                obj.BuildTree(currentNode.rightNodes);
            end
            
            
            
            
            
        end % build tree
        
        function h = visualize(obj)
            h = figure();
            axis equal;
            hold on;
            obj.visualizeNode(obj.rootNode);
            hold off;
            
            
            
        end
        
        function h = visualizeTree(obj)
           
            
            
            h = figure();            
            hold on;
            obj.visualizeTree_branch(obj.rootNode,0);
            hold off;
            
            end
        
            
        
        function  intersectTrees(obj, otherTree)
            % intersect two meshes from two trees
            obj.intersectTreeWithNode(otherTree.rootNode);
            obj.mesh{2} = otherTree.mesh{1};
            % Look for dupplicated nodes and remove them from the point
            % list
            
            npoints = size(obj.intersection_path,1);
            toRemove =[]; % two columns: left: the point to be removed, right the point to be replaced for in the conectivity
            for i=1:(npoints-1)
                for j=i+1:npoints
                    if (norm(obj.intersection_path(i,:)-obj.intersection_path(j,:))<10^-10)
                        % they are the same
                        toRemove = [toRemove ; j];
                        obj.newConectivity(obj.newConectivity==j)=i;
                    end
                end
            end
            % Remove duplicated points and update the conectivities
            toRemove = sort(unique(toRemove));
            while (numel(toRemove))
                obj.intersection_path(toRemove(1),:)=[];
                obj.newConectivity(obj.newConectivity>toRemove(1))=obj.newConectivity(obj.newConectivity>toRemove(1))-1;
                toRemove = toRemove-1;
                toRemove(1)=[];
            end
            
            % Remove single points
            obj.newConectivity(obj.newConectivity(:,1)==obj.newConectivity(:,2),:)=[];
            
        end
        
        
      
        function intersect_ray(obj,currentNode, points, vector)
            % This function deepens into the hierarchy of the current tree
            % until it finds the intersecting bounding box of the current
            % tree with the other node.
         
                
                if numel(currentNode.leftNodes)
                    obj.intersect_ray(currentNode.leftNodes, points, vector);
                end
                
                if numel(currentNode.rightNodes)
                    obj.intersect_ray(currentNode.rightNodes, points, vector);                    
                end
                
                if ~numel(currentNode.rightNodes) && ~numel(currentNode.rightNodes)
                   % newElements =[ currentNode.unique_triangles_id otherNode.unique_triangles_id ];
                    % if ( newElements(1)==newElements(2) || isequal(ismember(newElements, obj.computedIntersections),[1 1]))
                    %     return;
                    % end
                    % Compute intersection between two nodes (triangle vs
                    % triangle)
                   % intersection = currentNode.intersect_ray_with_triangle(points(points_that_can_be_inside>0,:), vector);
                    points_that_intersect_that_node = currentNode.intersect_ray_with_triangle(points, vector);
                    obj.computedIntersections(points_that_intersect_that_node) = obj.computedIntersections(points_that_intersect_that_node)+1;
                       

                end
            
        end
        
     
        
        function outMeshArray = generateOutputMesh(obj,booleanOp)
            % This function generated the resulting mesh. It can generate
            % three different meshes:
            %   booleanOp='difference'
            %   booleanOp='union'
            %   booleanOp='intersection'
            
            
            % For each triangle where intersection occurs, we have to
            % find the points on the triangle and build a constrained
            % triangularization
            
            %  PROCESS FOR MESH 1  ----------------------------------------
            % seed  is in the side I want to keep for this particular
            % boolean op.
            [outMesh1 seed1 ] = obj.processBoleanForOneMesh(1, booleanOp);
            if ~numel(seed1)
                outMeshArray=[];
                return;
            end
            [outMesh2 seed2 ] = obj.processBoleanForOneMesh(2, booleanOp);
            
            outMesh1.imposePositiveOrientation();
            outMesh2.imposePositiveOrientation();
            
            
            %%% ------------------------------ PROCESS MESHES -----------
            % for each mesh, give one label to the triangles at one point
            % of the intersection and a different labels to points at the
            % other side of the intersection.
            % This can be seen as a hole filling morphological operation.
            % To select a point in the positive side, we find the center of
            % the contour and project it on the surface positivewise
            
            outMesh1_segmenteda = meshFillHoles(outMesh1, seed1,'name','segmented');
            natr = outMesh1_segmenteda.find_attribute('segmented');
            
            npointstoremove1a = numel(nonzeros(outMesh1_segmenteda.attributes(natr).attribute_array));
            outMesh1_segmentedb = MeshType(outMesh1_segmenteda); % mesh copy
            
            outMesh1_segmentedb.attributes(natr).attribute_array= int8(...
                ~outMesh1_segmenteda.attributes(natr).attribute_array & ...
                ~outMesh1.attributes(outMesh1.find_attribute('roi')).attribute_array);
            npointstoremove1b = numel(nonzeros(outMesh1_segmentedb.attributes(natr).attribute_array));
            
            
            
            % remove the non-tagget vertices
            if (npointstoremove1a)
                outMesh1_segmenteda = outMesh1_segmenteda.removeVertex(find(outMesh1_segmenteda.attributes(natr).attribute_array));
            else
                % There is no point to remove, what should be removed is
                % triangles then. All triangle whose 3  triangles vertices
                % are tagged, in outMesh1
                % this will only happen either here or in the b tag . If it
                % happens in both, they are both the same (no intersection or total inclusion)
                taggedVertices = find(outMesh1.attributes(outMesh1.find_attribute('roi')).attribute_array);
                trianglesToRemove = [];
                for i=1:outMesh1.ntriangles
                    if ( ismember(outMesh1.triangles(i,1), taggedVertices) && ...
                            ismember(outMesh2.triangles(i,2), taggedVertices) && ...
                            ismember(outMesh1.triangles(i,3), taggedVertices))
                        trianglesToRemove = [trianglesToRemove; i];
                    end
                end
                
                outMesh1_segmenteda=outMesh1_segmenteda.removeTriangle(trianglesToRemove);
                
                
            end
            
            if (npointstoremove1b)
                outMesh1_segmentedb = outMesh1_segmentedb.removeVertex(find(outMesh1_segmentedb.attributes(natr).attribute_array));
            else
                % There is no point to remove, what should be removed is
                % triangles then. All triangle whose 3  triangles vertices
                % are tagged, in outMesh1
                % this will only happen either here or in the b tag . If it
                % happens in both, they are both the same (no intersection
                % or total inclusion)
                taggedVertices = find(outMesh1.attributes(outMesh1.find_attribute('roi')).attribute_array);
                trianglesToRemove = [];
                for i=1:outMesh1.ntriangles
                    if ( ismember(outMesh1.triangles(i,1), taggedVertices) && ...
                            ismember(outMesh1.triangles(i,2), taggedVertices) && ...
                            ismember(outMesh1.triangles(i,3), taggedVertices))
                        trianglesToRemove = [trianglesToRemove; i];
                    end
                end
                
                outMesh1_segmentedb=outMesh1_segmentedb.removeTriangle(trianglesToRemove);
            end
            
            
            % second mesh -----
            outMesh2_segmenteda = meshFillHoles(outMesh2, seed2,'name','segmented');
            natr = outMesh2_segmenteda.find_attribute('segmented');
            npointstoremove2a = numel(nonzeros(outMesh2_segmenteda.attributes(natr).attribute_array));
            
            outMesh2_segmentedb = MeshType(outMesh2_segmenteda); % mesh copy
            
            outMesh2_segmentedb.attributes(natr).attribute_array= int8(...
                ~outMesh2_segmenteda.attributes(natr).attribute_array & ...
                ~outMesh2.attributes(outMesh2.find_attribute('roi')).attribute_array);
            npointstoremove2b = numel(nonzeros(outMesh2_segmentedb.attributes(natr).attribute_array));
            
            % remove the non-tagget vertices
            if (npointstoremove2a)
                outMesh2_segmenteda = outMesh2_segmenteda.removeVertex(find(outMesh2_segmenteda.attributes(natr).attribute_array));
            else
                % There is no point to remove, what should be removed is
                % triangles then. All triangle whose 3  triangles vertices
                % are tagged, in outMesh1
                % this will only happen either here or in the b tag . If it
                % happens in both, they are both the same (no intersection
                % or total inclusion)
                taggedVertices = find(outMesh2.attributes(outMesh2.find_attribute('roi')).attribute_array);
                trianglesToRemove = [];
                for i=1:outMesh2.ntriangles
                    if ( ismember(outMesh2.triangles(i,1), taggedVertices) && ...
                            ismember(outMesh2.triangles(i,2), taggedVertices) && ...
                            ismember(outMesh2.triangles(i,3), taggedVertices))
                        trianglesToRemove = [trianglesToRemove; i];
                    end
                end
                
                outMesh2_segmenteda=outMesh2_segmenteda.removeTriangle(trianglesToRemove);
            end
            
            if (npointstoremove2b)
                outMesh2_segmentedb = outMesh2_segmentedb.removeVertex(find(outMesh2_segmentedb.attributes(natr).attribute_array));
            else
                % There is no point to remove, what should be removed is
                % triangles then. All triangle whose 3  triangles vertices
                % are tagged, in outMesh1
                % this will only happen either here or in the b tag . If it
                % happens in both, they are both the same (no intersection
                % or total inclusion)
                taggedVertices = find(outMesh2.attributes(outMesh2.find_attribute('roi')).attribute_array);
                trianglesToRemove = [];
                for i=1:outMesh2.ntriangles
                    if ( ismember(outMesh2.triangles(i,1), taggedVertices) && ...
                            ismember(outMesh2.triangles(i,2), taggedVertices) && ...
                            ismember(outMesh2.triangles(i,3), taggedVertices))
                        trianglesToRemove = [trianglesToRemove; i];
                    end
                end
                
                outMesh2_segmentedb=outMesh2_segmentedb.removeTriangle(trianglesToRemove);
                
            end
            
            
            
            
            % put together meshes??
            
            
            outMeshArray{1}=outMesh1_segmenteda;
            outMeshArray{2}=outMesh1_segmentedb;
            outMeshArray{3}=outMesh2_segmenteda;
            outMeshArray{4}=outMesh2_segmentedb;
        end
        
        
    end % end public functions
    
    
    methods(Access = private)
        
        
        
        
        function intersectTreeWithNode(obj, otherNode)
            % This function iterated and does "something" for all the nodes
            % of the other tree which lay inside the current tree
            
            if (obj.rootNode.intersect(otherNode))
                
                if (numel(otherNode.leftNodes))
                    obj.intersectTreeWithNode(otherNode.leftNodes)
                end
                
                if (numel(otherNode.rightNodes))
                    obj.intersectTreeWithNode(otherNode.rightNodes)
                end
                
                if (~numel(otherNode.rightNodes) && ~numel(otherNode.rightNodes))
                    % It is a leaf node: compute intersection
                    %intersection = obj.intersectNodeWithNode(obj.rootNode, otherNode);
                    obj.intersectNodeWithNode(obj.rootNode, otherNode);
                end
            end
            
        end
        
        function intersection = intersectNodeWithNode(obj, currentNode, otherNode)
            % This function deepens into the hierarchy of the current tree
            % until it finds the intersecting bounding box of the current
            % tree with the other node.
            intersection=NaN;
            
            if (currentNode.intersect(otherNode))
                
                if (numel(currentNode.leftNodes))
                    obj.intersectNodeWithNode(currentNode.leftNodes, otherNode);
                end
                
                if (numel(currentNode.rightNodes))
                    obj.intersectNodeWithNode(currentNode.rightNodes, otherNode);
                end
                
                if (~numel(currentNode.rightNodes) && ~numel(currentNode.rightNodes))
                    newElements =[ currentNode.unique_triangles_id otherNode.unique_triangles_id ];
                    % if ( newElements(1)==newElements(2) || isequal(ismember(newElements, obj.computedIntersections),[1 1]))
                    %     return;
                    % end
                    % Compute intersection between two nodes (triangle vs triangle)
                    intersection = currentNode.computeIntersection(otherNode);
                    if (~isnan(intersection))
                        obj.computedIntersections = [ obj.computedIntersections  ; newElements ];
                        %hold on; plot3(intersection(:,1),intersection(:,2), intersection(:,3),'r*-' ); hold off
                        obj.newConectivity = [obj.newConectivity  ; size(obj.intersection_path,1)+1 size(obj.intersection_path,1)+(size(intersection,1))];
                        obj.intersection_path = [obj.intersection_path; intersection];
                        
                    end
                end
            end
        end
        
        
        
        function visualizeNode(obj,currentNode)
            points = currentNode.bounds([1 3 5;
                1 3 6;
                1 4 5;
                1 4 6;
                2 3 5;
                2 3 6;
                2 4 5;
                2 4 6]);
            linerange = [1 2 4 3 1 5 6 2 6 8 4 8 7 3 7 5];
            mycolors = jet(obj.curr_max_depth+2);
            line(points(linerange,1),points(linerange,2),points(linerange,3),'color',mycolors(currentNode.depth+1,:));
            %,'LineWidth',obj.curr_max_depth -currentNode.depth...
            
            
            
            if (numel(currentNode.leftNodes))
                obj.visualizeNode(currentNode.leftNodes)
            end
            
            if (numel(currentNode.rightNodes))
                obj.visualizeNode(currentNode.rightNodes)
            end
            
            
        end
        
        function visualizeTree_branch(obj,currentNode,father_x)
            
            point_father = [father_x -currentNode.depth];
            
            if numel(currentNode.rightNodes) || numel(currentNode.leftNodes)
               plot(point_father(1),point_father(2),'o') ;
            end
            
            if (numel(currentNode.leftNodes))
                %point_son1 = [(currentNode.leftNodes.index_among_broders-0.5)/(currentNode.depth+1) -currentNode.leftNodes.depth];
                point_son1 = [(point_father(1))-1/(currentNode.depth+0.3) -currentNode.leftNodes.depth];
                
                line([point_father(1); point_son1(1)],[point_father(2); point_son1(2)]);
                
                obj.visualizeTree_branch(currentNode.leftNodes,point_son1(1));
                
           end
            
            if (numel(currentNode.rightNodes))
                point_son2 = [(point_father(1))+1/(currentNode.depth+0.3) -currentNode.rightNodes.depth];
                
                line([point_father(1); point_son2(1)],[point_father(2); point_son2(2)]);
                
                obj.visualizeTree_branch(currentNode.rightNodes,point_son2(1));
            end
            
            if ~numel(currentNode.rightNodes) && ~numel(currentNode.leftNodes)
               plot(point_father(1),point_father(2),'*') ;
            end
            
            
        end
        
        function [outMesh1 taggedVertex]= processBoleanForOneMesh(obj,meshIndex, booleanOp)
            % mesh index has to be 1 or 2
            % taggedVertex1 will be the one vertex on the positive side of
            % the mesh, which will be usable for filling holes later.
            
            
            indexOfTheOtherMesh = 1+mod(meshIndex,2); % if 2 then 1, if 1 then 2
            outMesh1 = MeshType(obj.mesh{meshIndex});
            
            taggedVertex = [];
            
            if ~numel(obj.computedIntersections)
                return;
            end
            
            intersectedTriangles1 = unique(obj.computedIntersections(:,meshIndex));
            for i=1:numel(intersectedTriangles1)
                % calculate local coordinate system
                corners = obj.mesh{meshIndex}.points( obj.mesh{meshIndex}.triangles(intersectedTriangles1(i),:),:);
                normal = cross(corners(2,:)-corners(1,:), corners(3,:)-corners(1,:))';
                normal = normal/norm(normal);
                
                [x_ y_] = vtkMathPerpendiculars(normal, pi/2);
                M = [x_ y_ normal];
                
                % Bring all points to the triangle plane and select
                % only the ones which are in plane
                mypoints3D = [obj.intersection_path ; corners ]';
                mypoints2D = M\mypoints3D ;
                
                % remove all points whose z component is far from the
                % tree original corners
                epsilon=10^-8;
                pointsInTriangle = abs((mypoints2D(3,:)-mypoints2D(3,end)))<epsilon;
                % Now, this can be not enough because if there is a
                % rectangular face, then the neighbor of the triangle
                % is in the same plane. We need to make sure it is IN
                % the triangle
                point_indices = find(pointsInTriangle(1:(end-3))); % remove 3 for the 3 corners
                corners2D = mypoints2D(1:2,(end-2):end);
                
                flag = PointInTriangle(mypoints2D(1:2,point_indices), corners2D(:,1),corners2D(:,2),corners2D(:,3));
                newPoints = point_indices(find(flag)); % of all the intersection points, these are in the triangle in the current mesh
                
                % see how this points are connected
                connectivity = obj.newConectivity(find(sum(ismember(obj.newConectivity,newPoints),2)==2),:);
                
                % sort the connectivity
                
                [ orderedNodes ordering] = obj.sortConnectivity(connectivity, newPoints);
                
                [coincidentPoints index_in_mesh index_in_orderedPoints]=intersect(outMesh1.points,obj.intersection_path(orderedNodes,:),'rows');
                clear coincidentPoints;
                
                
                if ~numel(taggedVertex)
                    
                    % find the normal to the triangle from the other
                    % surface, which will be used to determine edge
                    % orientation. If one facet has many segments, we only
                    % compute for one because the orientation of all the
                    % facet has to be consistent.
                    
                    % TODO:
                    % make that the orientation of the vector_of_facets is
                    % always the same. For that, the cross product has to
                    % be the same independently on if this is meshIndex=1
                    % or 2!!!!!!!!!!!!!1
                    
                    
                    triangle_intersected=[0 0]';
                    
                    triangle_intersected(meshIndex)=intersectedTriangles1(i);
                    connections = obj.computedIntersections(obj.computedIntersections(:,meshIndex)==intersectedTriangles1(i),:);
                    triangle_intersected(indexOfTheOtherMesh)= connections(1,indexOfTheOtherMesh);
                    % normal mesh 1
                    corners1 = obj.mesh{1}.points( obj.mesh{1}.triangles(triangle_intersected(1),:),:);                    
                    normal1 = cross(corners1(2,:)-corners1(1,:), corners1(3,:)-corners1(1,:))';
                    normal1 = normal1/norm(normal1);
                    % normal mesh 2
                    corners2 = obj.mesh{2}.points( obj.mesh{2}.triangles(triangle_intersected(2),:),:);
                    normal2 = cross(corners2(2,:)-corners2(1,:), corners2(3,:)-corners2(1,:))';
                    normal2 = normal2/norm(normal2);
                    
                    
                    %connections = obj.computedIntersections(obj.computedIntersections(:,meshIndex)==intersectedTriangles1(i),:);
                    %neighbour_triangle = connections(1,indexOfTheOtherMesh);
                    %cornersNeighbour = obj.mesh{indexOfTheOtherMesh}.points( obj.mesh{indexOfTheOtherMesh}.triangles(neighbour_triangle ,:),:);
                    %normal_neighbour = cross(cornersNeighbour(2,:)-cornersNeighbour(1,:), cornersNeighbour(3,:)-cornersNeighbour(1,:))';
                    %normal_neighbour = normal_neighbour/norm(normal_neighbour);
                    
                    %Now I have the path on the triangle. Let's decide the
                    %orientation I want. The orientation will be with respect
                    %to orderedNodes
                    
                    vector_of_facets =  cross(normal1,normal2);
                    intersection_pt = obj.intersection_path(newPoints(1),:);
                    
                    
                    
                    triangleCorners = obj.mesh{meshIndex}.triangles(intersectedTriangles1(i),:);
                    corners = obj.mesh{meshIndex}.points(triangleCorners,:); % points are rows
                    signedVectors = corners - ones(3,1)*intersection_pt;
                    sideOfSeed = sign(dot( ones(3,1)*vector_of_facets',signedVectors ,2));
                    positiveSeed = triangleCorners(sideOfSeed >0);
                    positiveSeed=positiveSeed(1);
                    negativeSeed = triangleCorners(sideOfSeed <0);
                    negativeSeed = negativeSeed(1);
                    
                    if (meshIndex==1)
                        if (strcmp(booleanOp,'difference'))
                            taggedVertex = positiveSeed;
                        elseif(strcmp(booleanOp,'intersection'))
                            taggedVertex = negativeSeed;
                        elseif (strcmp(booleanOp,'union'))
                            taggedVertex = positiveSeed;
                        elseif (strcmp(booleanOp,'clip'))
                            taggedVertex = negativeSeed;
                        else
                            disp('Error: not a good operator')
                            return;
                        end
                        
                        
                    else
                        
                        if (strcmp(booleanOp,'difference'))
                            taggedVertex = negativeSeed; % for second mesh
                        elseif(strcmp(booleanOp,'intersection'))
                            taggedVertex = negativeSeed; % for second mesh
                        elseif (strcmp(booleanOp,'union'))
                            taggedVertex = positiveSeed; % for second mesh
                        elseif (strcmp(booleanOp,'clip'))
                            taggedVertex = negativeSeed;
                        else
                            disp('Error: not a good operator')
                            return;
                        end
                        
                    end
                    
                end
                
                
                % redo delaunay
                mypoints2D = [ mypoints2D(1:2, orderedNodes) corners2D ] ;
                %make sure it is convex
                centre2D = mean(mypoints2D ,2);
                npts = numel(newPoints);
                ordered_segments = ordering([1 2]);
                for k=3:numel(ordering)
                    ordered_segments =[ordered_segments ; ordering([k-1 k])];
                end
                indices_of_points_to_be_added = (1:npts)+outMesh1.npoints;
                pointsOfTheNewTriangle = [ indices_of_points_to_be_added  obj.mesh{meshIndex}.triangles(intersectedTriangles1(i),:)];
                
                epsilon2=10^-12;
                
                % notsame = true;
                % while notsame
                %    niter = niter+1;
                mypoints2D_tmp = mypoints2D -centre2D*ones(1,size(mypoints2D,2));
                mypoints2D_tmp = mypoints2D_tmp*(1-epsilon2);
                mypoints2D(:,(npts+1):end) = mypoints2D_tmp(:,(npts+1):end)+centre2D*ones(1,3);
                clear DT;
                DT = DelaunayTri(mypoints2D(1,:)',mypoints2D(2,:)');
                %   triplot(DT);
                %    DT.Constraints=ordered_segments;
                %     notsame = max(DT.Triangulation(:))>numel(pointsOfTheNewTriangle);
                %     if niter>100
                %         disp('Error: Too many iterations in Delaunay!');
                %         return;
                %     end
                % end
                
                newTriangles =  pointsOfTheNewTriangle(DT.Triangulation);
                
                % in the new Triangles, replace the duplicated points, at
                % the positions
                % indices_of_points_to_be_added(index_in_orderedPoints) by
                % by index_in_mesh
                
                
                for j=1:numel(index_in_mesh)
                    newTriangles(newTriangles==indices_of_points_to_be_added(index_in_orderedPoints(j))) = index_in_mesh(j);
                    newTriangles(newTriangles>indices_of_points_to_be_added(index_in_orderedPoints(j)))= ...
                        newTriangles(newTriangles>indices_of_points_to_be_added(index_in_orderedPoints(j)))-1;
                    % remove the point and subtract 1 from the indices
                    % higher than that one
                    indices_of_points_to_be_added(indices_of_points_to_be_added>indices_of_points_to_be_added(index_in_orderedPoints(j)))=...
                        indices_of_points_to_be_added(indices_of_points_to_be_added>indices_of_points_to_be_added(index_in_orderedPoints(j)))-1;
                    
                end
                
                
                
                
                outMesh1.triangles(intersectedTriangles1(i),:)=newTriangles(1,:); % We do this way so that triangle indexing does notchange
                outMesh1.triangles = [outMesh1.triangles ; newTriangles(2:end,:)];
                outMesh1.ntriangles = size(outMesh1.triangles,1);
                
                points_to_be_added = obj.intersection_path(orderedNodes,:); % here we add newPoints instead of orderedPoints to be consistent with the loop above
                points_to_be_added(index_in_orderedPoints,:)=[];
                
                outMesh1.points = [outMesh1.points; points_to_be_added];
                outMesh1.npoints = size(outMesh1.points,1);
            end %
            
            
            
            % create an attribute which will have the points belonging to
            % the surface we are interested in
            at = AttributeType();
            at.attribute='field';
            at.name='roi';
            at.numComp=1; % between 1 and 4
            at.type='short';
            at.nelements=outMesh1.npoints; %number of tuples
            at.attribute_array=zeros(at.nelements,1);
            at.attribute_array((obj.mesh{meshIndex}.npoints+1):end)=1;
            
            ind = outMesh1.find_attribute(at.name);
            if (ind > 0)
                outMesh1.attributes(ind)=at;
            else
                n_attributes = numel(outMesh1.attributes);
                outMesh1.attributes(n_attributes+1)=at;
            end
        end
        
        
        
        
    end % end private functions
    
    methods(Static)
        function [ orderedNodes ordering] = sortConnectivity(connectivity, newPoints)
            firstNode=-1; % first node of the chain
            npts = numel(newPoints);
            for j=1:npts
                nocurrence(j)=numel(find(newPoints(j)==connectivity));
                if (nocurrence(j)==1)
                    firstNode=j;
                end
            end
            orderedNodes = newPoints(firstNode); % nodes ordered in the direction of the path but in an arbitrary sense
            ordering = firstNode;
            npts_=1;
            while (npts_<npts)
                [r c]= find(connectivity==orderedNodes(end));
                neighbor = connectivity(r,mod(c,2)+1);
                connectivity(r,:)=[];
                orderedNodes = [orderedNodes neighbor];
                if numel(neighbor)
                    ordering = [ordering find(newPoints==neighbor)]; % ordering so that newPoints(ordering ) = orderedNodes;
                else
                    % it can happen that in a triangular facet there
                    % are discontinuous lines. if ge are here it means
                    % that we have finished one segment and there are
                    % still more points to be connected separately
                    
                    % remove the points that have been already added
                    % and start again
                    newPoints= setdiff(newPoints,orderedNodes);
                    firstNode=-1; % first node of the chain
                    npts2 = numel(newPoints);
                    for j=1:npts2
                        nocurrence(j)=numel(find(newPoints(j)==connectivity));
                        if (nocurrence(j)==1)
                            firstNode=j;
                        end
                    end
                    orderedNodes = [orderedNodes newPoints(firstNode) ];
                    
                end
                npts_=npts_+1;
            end
        end
        
    end
    
    
end
