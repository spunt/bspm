function imageOut = meshClip(meshIn, meshStencil, varargin)
% Clips the input mesh, meshIn, using another mesh as a stencil
% 
%   Options:
%   none


dbg=false;
i=1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'thickness'))
        %n= varargin{i+2};
        %i = i+1;
    elseif(strcmp( varargin{i} , 'debug'))
        dbg= true;
    end
    i = i+1;
end

% build AABB tree

tree= AABB_tree(meshIn);
tree.BuildTree();

% input bounds

inputBounds = imageIn.GetBounds();
imageOut = ImageType(imageIn);

% For each leave of the tree, only if it its BB intersects with the BB of
% the input image, rasterize
imageOut = rasterize(tree.rootNode, imageOut,n);


end


function imageOut = rasterize(node, imageIn,n)
% This function iterated and does "something" for all the nodes
% of the other tree which lay inside the current tree

    imageOut = ImageType(imageIn);
    imageOut.data = imageIn.data;

bb = imageIn.GetBounds();
if (node.intersectBB(bb))
    % current node intersects with image, continue in this branch
    
    if (numel(node.leftNodes))
        % if there are branches on the left
        imageOut = rasterize(node.leftNodes, imageOut,n);
    end
    
    if (numel(node.rightNodes))
        % if there are branches on the right
        imageOut = rasterize(node.rightNodes, imageOut,n);
    end
    
    if (~numel(node.rightNodes) && ~numel(node.rightNodes))
        % It is a leaf node: compute intersection
        %intersection = obj.intersectNodeWithNode(obj.rootNode, otherNode);
        imageOut = rasterizeNode(node, imageOut,n);
    end
end

end

function imageOut = rasterizeNode(currentNode, imageIn,n)
    % rasterize the current node into the image
    
    imageOut = ImageType(imageIn);
    imageOut.data = imageIn.data;
       
    ntriangles = currentNode.mesh.ntriangles;
    if (ntriangles > 1)
       disp('ERROR: at the end of the branch there should be only one triangle!') 
       return ;
    end
    
    % Find coordinate system of the triangle
    
    normal = currentNode.mesh.GetTriangleNormal(1);
    [x y ] = vtkMathPerpendiculars(normal,pi/2);
    
    M = [x y normal(:)];
    
    
    % convert imagePoints to the trangle plane
    [x y z] = ndgrid(1:imageIn.size(1), 1:imageIn.size(2), 1:imageIn.size(3));
    imagePoints = imageIn.GetPosition([x(:) y(:) z(:)]');
    
    imagePoints_transformed = M \ imagePoints;
    
    % see if the points are in the triangle
    triangleCorners = currentNode.mesh.points(currentNode.mesh.triangles(1,:),:)';
    triangleCorners_transformed = M \ triangleCorners;
    
    flag = PointInTriangle(imagePoints_transformed(1:2,:) , triangleCorners_transformed(1:2,1),triangleCorners_transformed(1:2,2),triangleCorners_transformed(1:2,3));
    
    %pointsCloser = abs(imagePoints_transformed(3,:)-triangleCorners_transformed(3,1))<imageIn.spacing(3);
    pointsCloser = abs(imagePoints_transformed(3,:)-triangleCorners_transformed(3,1))<n*norm(imageIn.spacing);
    
    pointsInTriangle = flag & pointsCloser;
    pointsInTriangle = reshape(pointsInTriangle, imageIn.size' );
    
    imageOut.data(pointsInTriangle)=1;

end