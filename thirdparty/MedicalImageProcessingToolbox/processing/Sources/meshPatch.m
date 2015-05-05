function [mesh, addedPoints, contour, points_spherical_all] = meshPatch( patch,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% A patch has the following members:
%    patch.circles: [1x18 struct]
%    patch.frustum: [1x1 struct]
%    patch.radius: 68.9241

mesh = [];
angle_res = 45;
points = [];
N = 22;
 for i=1:size(varargin,2)
    if (strcmp(varargin{i},'resolution'))
        N=varargin{i+1};
    end                
end



for j=1:numel(patch.circles)
    % for each circle in one patch
    m_patches  = circle3DMesh(patch.circles(j).centre,patch.circles(j).radius,...
        patch.circles(j).direction,'resolution',angle_res,'theta',patch.circles(j).theta);
    if j>1
        n(1) = norm(points(:,end)-m_patches.points(:,1));
        n(2) = norm(points(:,end)-m_patches.points(:,end));
        n(3) = norm(points(:,1)-m_patches.points(:,1));
        n(4) = norm(points(:,1)-m_patches.points(:,end));
        
        [~,idx]=min(n);
        switch idx
            case 1
                points =  [points  m_patches.points];
            case 2
                points =  [points  m_patches.points(:,end:-1:1)];
            case 3
                points =  [m_patches.points(:,end:-1:1) points  ];
            case 4
                points =  [m_patches.points points  ];
            otherwise
                disp('ERROR: bad value');
        end
    else
        points = [points  m_patches.points];        
    end
end

contour.cartesian = points;

% convert the points to spherical 
M = [ patch.frustum.directions patch.frustum.centre; 0 0 0 1];
points_centred = M\[points; ones(1,size(points,2))];
points_spherical = pointFromCartesianToSpherical(points_centred(1:3,:));
contour.spherical = points_spherical;
% make an image out of the contour

sp = max((max(points_spherical(1:2,:)')-min(points_spherical(1:2,:)'))/N);
N = ceil((max(points_spherical(1:2,:)')-min(points_spherical(1:2,:)'))/sp)+1+2; % 2 for borders
im = zeros(N);
origin = min(points_spherical(1:2,:)');
zvalue = mean(points_spherical(3,:));
ix = round((points_spherical(1,:)-origin(1))/sp)+2;
iy = round((points_spherical(2,:)-origin(2))/sp)+2;
ind = sub2ind(N,ix,iy);
im(ind)=1;

im = imfill(im);
er1=1;
[x,y,z] = ndgrid(-er1:er1,-er1:er1,-er1:er1);
se = (x/er1).^2 + (y/er1).^2 + (z/er1).^2 <= 1;
im = imerode(im,se);
% now convert nonzero values to 2d coordinates
ind = find(im);
if ~numel(ind)
   mesh = MeshType();
   addedPoints=NaN;
   contour=NaN;
   points_spherical_all=NaN;
   return;
end
[ix,iy] = ind2sub(N,ind);
points_spherical2(1,:) = (ix-2)*sp+origin(1);
points_spherical2(2,:) = (iy-2)*sp+origin(2);

% calculate triangularization
points_tr = [points_spherical(1:2,:) points_spherical2(1:2,:)];
% constraint border
c(:,1) = (1:size(points_spherical,2))';
c(:,2) = [2:size(points_spherical,2) 1]';
dt = DelaunayTri(points_tr(1,:)',points_tr(2,:)',c);
io = dt.inOutStatus();

triangles = dt.Triangulation(io,:);
points_spherical_all(1:2,:)=dt.X';
points_spherical_all(3,:) = zvalue*ones(1,size(dt.X,1));
%  convert to world coordinates
points_centred_all = pointFromSphericalToCartesian(points_spherical_all);
points_all = M * [points_centred_all; ones(1,size(points_centred_all,2))];

addedPoints =  M *[pointFromSphericalToCartesian([ points_spherical2(1:2,:) ; zvalue*ones(1,size(points_spherical2,2))]) ; ones(1,size(points_spherical2,2)) ];

mesh = MeshType();
mesh.triangles = triangles;
mesh.points = points_all(1:3,:)';
mesh.ntriangles = size(mesh.triangles,1);
mesh.npoints = size(mesh.points,1);


end

