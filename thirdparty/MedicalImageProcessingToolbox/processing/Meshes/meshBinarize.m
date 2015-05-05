function imageOut = meshBinarize(meshIn, imageIn, varargin)
% converts a mesh into a binary image (rasterization), on the grid defined
% 
%   Options:
%   'thickness' , n     ->  n is the number of voxel sizes at each size of
%   the mesh that will be painted (by default, n=1/2)

n=1/2;
dbg=false;
i=1;
enlarge=0;
MAX_CHUNK_SIZE = 50;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'thickness'))
        n= varargin{i+1};
        i = i+1;
    elseif (strcmp( varargin{i} , 'enlarge'))
        enlarge= varargin{i+1};
        i = i+1;
    elseif(strcmp( varargin{i} , 'debug'))
        dbg= true;
    elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end
    i = i+1;
end

% build AABB tree

meshOut = MeshType(meshIn);
meshOut.triangles =meshIn.triangles;
imageOut = ImageType(imageIn);

if enlarge
    centerpoint = mean(meshIn.points);
    if enlarge <0
        % this only works for planes! Compute intersection of bounding box
        % and plane
        plane.normal = meshIn.GetNormalAtFaces(1);
        plane.origin = meshIn.points(meshIn.triangles(1),:);
        % transform positions to align with normal
        [x y] = vtkMathPerpendiculars(plane.normal',pi/2);
        M = [x y plane.normal(:); 0 0 0];
        M = [M [plane.origin(:) ; 1] ];
        NCHUNKS = ceil(imageIn.size/MAX_CHUNK_SIZE);
        chunked_size = ceil(imageIn.size./NCHUNKS)';
    
        [ix iy iz]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1);
        intervals = [ix(:) iy(:) iz(:)];
        clear ix iy iz;
        for i=1:size(intervals,1)
            ranges([1 3 5]) = intervals(i,:).*chunked_size+1;
            ranges([2 4 6]) = min([(intervals(i,:)+[1 1 1]).*chunked_size ; imageIn.size']);
        
            ranges_size = ranges([2 4 6])-ranges([1 3 5])+[1 1 1];
            % generate all the indexes of the target image
            [x y z] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6));
            positions = imageIn.GetPosition([x(:) y(:) z(:)]');
            clear x y z;
            
            tx_positions = M \ [positions ; ones(1,size(positions,2)) ];
            
            imageOut.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(abs(tx_positions(3,:)')<n,ranges_size);
        end
        return;
        
    else
        for i=1:meshIn.npoints
            meshOut.points(i,:)=(meshIn.points(i,:)-centerpoint)*enlarge+ centerpoint;
        end
    end
else
    meshOut.points=meshIn.points;
end

tree= AABB_tree(meshOut);
tree.BuildTree();

% input bounds




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
    
    if (currentNode.mesh.npoints < 3)
        disp('ERROR: This triangle has less than 3 points!') 
        return ;
    end
        
    
    % Find coordinate system of the triangle
    
    normal = currentNode.mesh.GetTriangleNormal(1);
    [x y ] = vtkMathPerpendiculars(normal(:),pi/2);
    
    M = [x y normal(:)];
    
    
    % convert imagePoints to the trangle plane
    [x y z] = ndgrid(1:imageIn.size(1), 1:imageIn.size(2), 1:imageIn.size(3));
    imagePoints = imageIn.GetPosition([x(:) y(:) z(:)]');
    
    imagePoints_transformed = M \ imagePoints;
    
    % see if the points are in the triangle
    triangleCorners = currentNode.mesh.points(currentNode.mesh.triangles(1,:),:)';
    triangleCorners_transformed = M \ triangleCorners;
    
    flag = PointInTriangle(imagePoints_transformed(1:2,:) , triangleCorners_transformed(1:2,1),triangleCorners_transformed(1:2,2),triangleCorners_transformed(1:2,3),'maxChunkSize',150);
    
    %pointsCloser = abs(imagePoints_transformed(3,:)-triangleCorners_transformed(3,1))<imageIn.spacing(3);
    pointsCloser = abs(imagePoints_transformed(3,:)-triangleCorners_transformed(3,1))<n*norm(imageIn.spacing);
    
    pointsInTriangle = flag & pointsCloser;
    pointsInTriangle = reshape(pointsInTriangle, imageIn.size' );
    
    imageOut.data(pointsInTriangle)=1;

end