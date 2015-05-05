function out = prismMesh(roi_file, varargin)
% generate a prism  from a roi, using the roi as base section and extruding
% along the normal to the roi plane. Default height is 100  mm

height = 100; %mm
i=1;
closeLids = false;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'height'))
        height= varargin{i+1};
        i = i+1;
    elseif(strcmp( varargin{i} , 'closeLids'))
        closeLids= true;
    elseif(strcmp( varargin{i} , 'debug'))
        dbg= true;
    end
    i = i+1;
end

%% parameters

%% Read roi 
% ------------------------------------------- read roi

[roinormal, roiorigin, points2D, points] =  read_roi(roi_file);

points = points(1:3,:);

points_above = points + height*roinormal*ones(1,size(points,2));

outpoints = [points points_above ]';

outpoints = outpoints - (height/2*roinormal*ones(1,size(outpoints,1)))';
%out.points = [points points_above points_below]';

% build topology

npoints = size(points,2);

% make the triangles
outtriangles = [];

for i=1:npoints-1
    
    % above
    outtriangles =[outtriangles; i i+1 i+npoints; i+1 i+npoints+1 i+npoints ];
    % below
    
        
end


out = MeshType(size(outpoints,1),size(outtriangles,1)); 

out.points = outpoints;
out.triangles = outtriangles;

if closeLids 
    warning off;
    c = [1:size(points2D,2) ; 2:size(points2D,2) 1]';
    dt = DelaunayTri(points2D(1,:)',points2D(2,:)',c);
    io = dt.inOutStatus();
    lidtriangles = dt.Triangulation(io,:);
    lidtriangles2 = lidtriangles +npoints;    
    out.triangles =[out.triangles ;lidtriangles;lidtriangles2];
    warning on;
end

out.npoints = size(outpoints,1);
out.ntriangles= size( out.triangles,1);

end



