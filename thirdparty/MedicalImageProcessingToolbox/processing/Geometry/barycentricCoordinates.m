function barycoord= barycentricCoordinates(p, a,b,c,varargin)
% returns the barycentriuc coordinates of p. The user must guarantee that p
% is in the triangle
    
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
    
    barycoord = zeros(3,size(p,2));

    
    if norm(a-c)<1E-10 || norm(a-b)<1E-10
        return;
    end
    
    for i=1:numel(intervals)
        range(1) = intervals(i)*chunked_size+1;
        range(2) = min([(intervals(i)+1)*chunked_size  size(p,2)]);
        T = [a(1)-c(1) b(1)-c(1); a(2)-c(2) b(2)-c(2)];
        lambdas = T\(p-c*ones(1,size(p,2)));
        l3 = ones(1,size(p,2))-sum(lambdas,1);
        %barycoord(:,range(1):range(2)) = [a b c]\p;
        barycoord(:,range(1):range(2)) = [lambdas; l3];
    end
        
end
