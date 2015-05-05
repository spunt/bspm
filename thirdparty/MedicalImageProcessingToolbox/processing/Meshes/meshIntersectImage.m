function flag = meshIntersectImage(meshIn, imageIn, varargin)
% returns true if the mesh and the image intersect at voxels where the
% image is nonzero

n=1/2;
dbg=false;
i=1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'debug'))
        dbg= true;
    elseif (strcmp( varargin{i} , 'thickness'))
        n = varargin{i+1};
        i=i+1;
    end
    i = i+1;
end

% build AABB tree

tree= AABB_tree(meshIn);
tree.BuildTree();

% input bounds



% For each leave of the tree, only if it its BB intersects with the BB of
% the input image, rasterize
flag = find_intersection(tree.rootNode, imageIn,n);


end


function flag = find_intersection(node, imageIn,n)
% This function iterated and does "something" for all the nodes
% of the other tree which lay inside the current tree
flag=0;
bb = imageIn.GetBounds();
bb=bb(:);
if (node.intersectBB(bb))
    % current node intersects with image, continue in this branch
    
    if (numel(node.leftNodes))
        % if there are branches on the left
        flag = find_intersection(node.leftNodes, imageIn,n);
        if flag
            return;
        end
    end
    
    if (numel(node.rightNodes))
        % if there are branches on the right
        flag = find_intersection(node.rightNodes, imageIn,n);
        if flag
            return;
        end
    end
    
    if (~numel(node.rightNodes) && ~numel(node.rightNodes))
        % It is a leaf node: compute intersection
        %intersection = obj.intersectNodeWithNode(obj.rootNode, otherNode);
        flag = find_intersection_node(node, imageIn,n);
        
    end
  
end

end

function flag = find_intersection_node(currentNode, imageIn,n)
      flag = false;
        
    ntriangles = currentNode.mesh.ntriangles;
    if (ntriangles > 1)
       disp('ERROR: at the end of the branch there should be only one triangle!') 
       return ;
    end
    
    if (currentNode.mesh.npoints < 3)
        disp('ERROR: This triangle has less than 3 points!') 
        return ;
    end
        
    
    % Find coordinate system of the triangle
    
    normal = currentNode.mesh.GetTriangleNormal(1);
    [x y ] = vtkMathPerpendiculars(normal(:),pi/2);
    
    M = [x y normal(:)];
    
    
    % convert nonzero imagePoints to the trangle plane
    nonzero_indices = find(abs(imageIn.data)>0);
    if ~numel(nonzero_indices)
        return ;
    end
    [x y z] = ind2sub(imageIn.size',nonzero_indices);
    imagePoints = imageIn.GetPosition([x(:) y(:) z(:)]');    
    imagePoints_transformed = M \ imagePoints;
    
    % see if the points are in the triangle
    triangleCorners = currentNode.mesh.points(currentNode.mesh.triangles(1,:),:)';
    triangleCorners_transformed = M \ triangleCorners;
    
    fl = PointInTriangle(imagePoints_transformed(1:2,:) , triangleCorners_transformed(1:2,1),triangleCorners_transformed(1:2,2),triangleCorners_transformed(1:2,3),'maxChunkSize',150);
    
    %pointsCloser = abs(imagePoints_transformed(3,:)-triangleCorners_transformed(3,1))<imageIn.spacing(3);
    pointsCloser = abs(imagePoints_transformed(3,:)-triangleCorners_transformed(3,1))<n*norm(imageIn.spacing);
    pointsInTriangle = fl & pointsCloser;
    flag = numel(nonzeros(pointsInTriangle))>0;
    
    
    

end