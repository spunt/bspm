% test_generateCone
function out = prismImage(roi_file, ref_im,varargin)
% generate a prism  from a roi, using the roi as base section and extruding
% along the normal to the roi plane

%% parameters
MAX_CHUNK_SIZE = 50;
MAX_CHUNK_SIZE = 50;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'debug'))
        dbg = true;
    elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end
end
out = ImageType(ref_im);
th = sqrt(sum(ref_im.spacing.^2))/2;
npieces=1;


%% Read roi
% ------------------------------------------- read roi

[roinormal roiorigin points] =  read_roi(roi_file);

roi.normal = roinormal;
roi.origin = roiorigin;

[x y]=vtkMathPerpendiculars(roi.normal,pi/2);

M=eye(4);
M(1:3,1:3) = [x(:) y(:) roi.normal(:)];
M(1:3,4) = roi.origin;


npoints = size(points,2);

%% The image will be filled in

NCHUNKS = ceil(ref_im.size/MAX_CHUNK_SIZE);


chunked_size = ceil(ref_im.size./NCHUNKS)';

[ix iy iz]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1);
intervals = [ix(:) iy(:) iz(:)];
clear ix iy iz;
for i=1:size(intervals,1)
    ranges([1 3 5]) = intervals(i,:).*chunked_size+1;
    ranges([2 4 6]) = min([(intervals(i,:)+[1 1 1]).*chunked_size ; ref_im.size']);
    
    ranges_size = ranges([2 4 6])-ranges([1 3 5])+[1 1 1];
    % generate all the indexes of the target image
    [x y z] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6));
    
    ref_positions = ref_im.GetPosition([x(:) y(:) z(:)]');
    clear x y z;
    ref_positions_2D = M \ [ref_positions; ones(1,size(ref_positions,2))];
    clear ref_positions;
    % The origin of the system is the centroid of the points
    centroid = mean(points,2);
    % make the triangles
    
    for j=1:npoints-1
        % triangle i is defined by point i, point i+1 and centroid. Seee if the
        % ref positions are inside it
        A = points(1:2,j);
        B = points(1:2,j+1);
        
        pointsInside = PointInTriangle(ref_positions_2D(1:2,:), A,B,centroid(1:2));
        pointsInside = reshape(pointsInside,ranges_size);
        out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = ...
            double(out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) | pointsInside>0);
        
    end
    clear ref_positions_2D  pointsInside ;
end

end



