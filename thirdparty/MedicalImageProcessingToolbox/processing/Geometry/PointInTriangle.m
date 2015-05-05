function flag = PointInTriangle(p, a,b,c,varargin)
% returns true if the point p is inside the triangle defined by a,b,c (in 2D)
% p, a b and c are column points
    
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

    
    if norm(a-c)<1E-10 || norm(a-b)<1E-10
        return;
    end
    
    for i=1:numel(intervals)
        range(1) = intervals(i)*chunked_size+1;
        range(2) = min([(intervals(i)+1)*chunked_size  size(p,2)]);
        
        flag(range(1):range(2)) = SameSide(p(:,range(1):range(2)),a, b,c) .* SameSide(p(:,range(1):range(2)),b, a,c) .* SameSide(p(:,range(1):range(2)),c, a,b);
    end
        
end

function flag = SameSide(p1,p2, a,b)
% This function and the next one are taken from http://www.blackpawn.com/texts/pointinpoly/default.html 

    npts = size(p1,2);
    p1 = [p1; zeros(1,npts)];
    %p2Large = p2*ones(1,npts);
    
    aLarge = [a; 0]*ones(1,npts);
    bLarge = [b; 0]*ones(1,npts);
    
    cp1 = cross(bLarge-aLarge, p1-aLarge,1);
    cp2 = cross([b-a; 0], [p2; 0]-[a; 0],1);
    cp2Large = cp2*ones(1,npts);
    
    epsilon = 10^-10;
    
    flag = dot(cp1,cp2Large,1) >= (0-epsilon);
        
end


