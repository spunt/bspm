function [meshOut, borders] = closeMeshWithPlanes(meshIn_, varargin)

% only works for two planes

labelIndex = [];
th=0.001; % a percent
npts = 5; % min number of points in a border
i=1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'th'))
        th= varargin{i+1};
        i = i+1;
    elseif (strcmp( varargin{i} , 'n'))
        npts= varargin{i+1};
        i = i+1;
    end
    i = i+1;
end
meshIn = MeshType(meshIn_);
meshIn.detectBorders();
border_indices = find(meshIn.attributes(2).attribute_array);
% viewMesh(meshIn); hold on; plotpoints(meshIn.points(border_indices,:),'*'); hold off


% order borderpoints:

ordered_indices =border_indices(1);

while numel(ordered_indices)<numel(border_indices)
    %distance to last added point
    candidates = setdiff(border_indices, ordered_indices);
    
    dist = sum((meshIn.points(candidates,:)-ones(numel(candidates),1)*meshIn.points(ordered_indices(end),:)).^2,2);
    [~,new_index]=min(dist);
    ordered_indices = [ordered_indices candidates(new_index)];    
end
borderpoints = meshIn.points(ordered_indices,:);


indices = [ 1:numel(ordered_indices) 1:numel(ordered_indices) 1:numel(ordered_indices) 1:numel(ordered_indices) 1:numel(ordered_indices)];

first=numel(ordered_indices)*2;
last = numel(ordered_indices)*3;

error=Inf;
while error >th
    
    idx = first:last;
    
    
    
    centroid = mean(borderpoints(indices(idx),: ));
    
    S=borderpoints(indices(idx),: )-ones(size(borderpoints(indices(idx),: ),1),1)*centroid;
    [D,V]=eig(S'*S);
    
    plane.normal = D(:,1);
    M = D(:,[2 3 1]);
    
    tx_points = M\S';
    
   gaussian_radius = sqrt(std(tx_points(1,:)).^2+ std(tx_points(1,:)).^2);
    error = max(abs(tx_points(3,:)))/gaussian_radius;
    
    last=last-1;
    
end

while error <=th
    
    idx = first:last;
    
    
    
    centroid = mean(borderpoints(indices(idx),: ));
    
    S=borderpoints(indices(idx),: )-ones(size(borderpoints(indices(idx),: ),1),1)*centroid;
    [D,V]=eig(S'*S);
    
    plane.normal = D(:,1);
    M = D(:,[2 3 1]);
    
    tx_points = M\S';
    gaussian_radius = sqrt(std(tx_points(1,:)).^2+ std(tx_points(1,:)).^2);
    error = max(abs(tx_points(3,:)))/gaussian_radius;
    
    
    first=first-1;
    
    
end


border1 = indices(idx(2:end));
nremaining = numel(ordered_indices)-numel(border1);
border2 = indices(idx(end):idx(end)+nremaining+1 );

% calculate centroid of the line connecting both borders

centroid = mean(borderpoints( border1([1 end]),:));

borders = ordered_indices(border1([1 end]));

% new topology

meshOut =MeshType(meshIn);
meshOut.points = [meshIn.points; centroid];
meshOut.npoints = size(meshOut.points,1);

t=[];

for i=1:numel(border1)-1
    t = [t; ordered_indices([border1(i) border1(i+1)])  meshOut.npoints];
end

for i=1:numel(border2)-1
    t = [t; ordered_indices([border2(i) border2(i+1)]) meshOut.npoints];
end

meshOut.triangles= [meshIn.triangles; t];
meshOut.ntriangles = size(meshOut.triangles,1);




