function flag = PointInFrustum(p, frustum,varargin)
% returns true if the point p is inside the frustum (in 3D)
% p are column points
    
MAX_CHUNK_SIZE = 50;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end    
end



    NCHUNKS = ceil(size(p,2)/MAX_CHUNK_SIZE.^3);
    chunked_size = ceil(size(p,2)./NCHUNKS);
    intervals = 0:(NCHUNKS-1);
    
    flag = zeros(1,size(p,2));
    
    for i=1:numel(intervals)
        range(1) = intervals(i)*chunked_size+1;
        range(2) = min([(intervals(i)+1)*chunked_size  size(p,2)]);
        
        flag(range(1):range(2)) = isInside(p(:,range(1):range(2)),frustum) ;
    end
        
end

function flag = isInside(p,frustum)

    npts = size(p,2);
    flag = ones(1,npts);
     M = [frustum.directions frustum.centre(:) ; 0 0 0 1];
    points_in_frustum_coordinates = M\p;
    
    radiuses = norm(points_in_frustum_coordinates(1:3,:));
    flag(radiuses<frustum.radiuses(1))=0;
    flag(radiuses>frustum.radiuses(2))=0;
    
    theta_frustum = [atan2(points_in_frustum_coordinates(1,:),points_in_frustum_coordinates(3,:))
                    atan2(points_in_frustum_coordinates(2,:),points_in_frustum_coordinates(3,:))]*180/pi; % first row: angle1, second row: angle2
    
    a1 = angleInRange(theta_frustum(1,:),frustum.angles([1 2]),1E-03);
    a2 = angleInRange(theta_frustum(2,:),frustum.angles([3 4]),1E-03);
    inAngle = intersect(a1,a2);
    outAngle = setdiff(1:npts,inAngle);
    flag(outAngle)=0;
        
end


