function m  = circle3DMesh(c,r,d,varargin)
% m = circle3DMesh(c,r,d)
% m = circle3DMesh(c,r,d,options)
%
% creates a MeshType object of a sphere of radius r and center c, within a
% plane normal to vector d(:,3). d(:,1:2) are the other two axes, and theta
% is measured with respect to axis d(:,1)
%
% options       meaning     default
% -------       -------     -------
%
% 'resolution'      sphereRes   16
% 'theta'           angle       [-180 180]



resolution = 8; % half resolution for phi than for theta
theta = [0 360].*pi/180;
nOptions = size(varargin,2);
% Argument reading

if (nOptions > 0)
    if (rem(nOptions,2)~=0)
        disp('Error in the arguments, please check the option list');
        return;
    end
    i=1;
    while(i<=nOptions)
        
        if (strcmp(varargin{i},'resolution'))
            resolution = varargin{i+1};
            i = i+2;
        elseif  (strcmp(varargin{i},'theta'))
            theta = varargin{i+1}*pi/180;
            i = i+2;
        end
    end
end

%theta(1)=max(theta(1),0);
%theta(2)=min(theta(2),360);

closed = true;
if (theta(2)-2*pi< theta(1))
    closed = false;
end

% --------------------------


thetaResolution = resolution;


numPts =  2*resolution+2;
numPolys = 2*resolution;



% Create circle

xv = d(:,1);
yv = d(:,2);
n = d(:,3);

if theta(1)<theta(2)
    theta = theta(1):((theta(2)-theta(1))/(resolution)):theta(2);
    if (size(theta,2)==resolution)
        % add last point
        theta(end+1)=theta(end)+((theta(2)-theta(1))/(resolution));
    end
else
    theta = theta(1):((theta(2)-theta(1)+2*pi)/(resolution)):(theta(2)+2*pi);
    if (size(theta,2)==resolution)
        % add last point
        theta(end+1)=theta(end)+((theta(2)-theta(1)+2*pi)/(resolution));
    end
end



% points have to go from x to y

points = c*ones(1,resolution+1) + ...
    r*([1 1 1]'*cos(theta)).*(xv*ones(1,resolution+1)) + ...
    r*([1 1 1]'*sin(theta)).*(yv*ones(1,resolution+1));


polygons = [];

for s=1:resolution
    polygons = [polygons
        s s+1 s+1];
end


% Create conectivity



% create mesh
m = MeshType(size(points,1),size(polygons,1));
m.points =points;
m.triangles = polygons;


% just for testing

end
