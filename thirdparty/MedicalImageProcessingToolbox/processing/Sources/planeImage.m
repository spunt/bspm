
function out = planeImage(p, n, ref_im,varargin)
%lets say that the plane is defined by the point p and the
%normal vector n.

skeletonize_plane=false;
MAX_CHUNK_SIZE = 50;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'debug'))
        dbg = true;
    elseif (strcmp(varargin{i},'skeletonize'))
        skeletonize_plane= true;
        i=i+1;
    elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end
    
end

%% parameters

out = ImageType(ref_im);
th = sqrt(sum(ref_im.spacing.^2))/2;


NCHUNKS = ceil(out.size/MAX_CHUNK_SIZE);
chunked_size = ceil(out.size./NCHUNKS)';

[ix iy iz]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1);
intervals = [ix(:) iy(:) iz(:)];
clear ix iy iz;
for i=1:size(intervals,1)
    ranges([1 3 5]) = intervals(i,:).*chunked_size+1;
    ranges([2 4 6]) = min([(intervals(i,:)+[1 1 1]).*chunked_size ; out.size']);
    
    ranges_size = ranges([2 4 6])-ranges([1 3 5])+[1 1 1];
    % generate all the indexes of the target image
    [x y z] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6));
    
    %[ x y z] = ndgrid(1:out.size(1), 1:out.size(2), 1:out.size(3));
    positions = out.GetPosition([x(:) y(:) z(:)]');
    clear x y z;
    planeToPoint =  p*ones(1,size(positions,2)) - positions;
    clear  positions;
    distancesToPlane = abs(n'*planeToPoint);
    clear planeToPoint;
    distance_th = sqrt(sum(ref_im.spacing.^2))/2; % world coordinates
    outputValues = double(distancesToPlane < distance_th  );
    clear distancesToPlane;
    %out.data(distancesToPlane<distance_th)=1;
    out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(outputValues,ranges_size);
    clear outputValues;
end

bounds = out.GetBounds(0.9,true);
out = cropImage(out,bounds);
if (skeletonize_plane)
    out.data = skeletonize3D(out.data);
end


end
