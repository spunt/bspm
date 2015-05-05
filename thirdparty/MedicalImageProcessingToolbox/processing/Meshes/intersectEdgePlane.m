function [point indices] = intersectEdgePlane(edge, plane, varargin)
%INTERSECTEDGEPLANE Return intersection point between a plane and a edge
%
%   PT = intersectEdgePlane(edge, PLANE) return the intersection point of
%   the given edge and the given plane.
%   PLANE : plane.point [x0 y0 z0]', plane.normal [nx ny nz]'
%   edge :  .p1 and .p2
%   PT :    [xi yi zi]'
%   If EDGE and PLANE are parallel, return [NaN NaN NaN]'.
%   If EDGE.px (or PLANE) is a matrix with 3 rows and N columns, result
%   is an array of points with N cols and N rows.
%
%
%
%
%   ---------
%   Alberto Gomez, inspired by
%   author : David Legland
%   INRA - TPV URPOI - BIA IMASTE
%   created the 24/04/2007 from intersectLinePlane.
%
%
%   HISTORY
%
%   14/05/2012 Modified by Alberto Gomez
%   17/06/2011 E. J. Payton - fixed indexing error that caused incorrect
%              points to be returned

tol = 1e-14;
if ~isempty(varargin)
    tol = varargin{1};
end

np = size(plane.point, 2);
ne = size(edge.p1, 2);

% unify sizes of data
if np ~= ne
    if ne == 1;
        % one edge and many planes
        edge.p1 = edge.p1(:,ones(1,np));
        edge.p2 = edge.p2(:,ones(1,np));
    elseif np == 1
        % one plane possible many edges
        plane.point = plane.point(:,ones(1, ne));
        plane.normal = plane.normal(:,ones(1, ne));
    else
        % N planes and M edges, not allowed for now.
        error('Should have the same number of planes and edges');
    end
end

% initialize empty arrays
point = zeros(3, size(plane.point, 2));
t = zeros(1, size(plane.normal,2)); % line parameter



% create line supporting edge
line.directorVector = edge.p2 - edge.p1;
line.length =  sqrt(line.directorVector(1,:).^2+line.directorVector(2,:).^2+line.directorVector(3,:).^2);
%for i=1:3
%    line.directorVector(i,:) = line.directorVector(i,:)./line.length;
%end

% get indices of edge and plane which are parallel
par = abs(dot(plane.normal, line.directorVector, 1)) < tol;
point(:,par) = NaN;
t(par) = NaN;

% difference between origins of plane and edge
dp = plane.point - edge.p1;

% relative position of intersection on line
%t = dot(n(~par,:), dp(~par,:), 2)./dot(n(~par,:), line(~par,4:6), 2);
t(~par) = dot(plane.normal(:,~par), dp(:,~par), 1) ./ dot(plane.normal(:,~par), line.directorVector(:,~par), 1);

% compute coord of intersection point
%point(~par, :) = line(~par,1:3) + repmat(t,1,3).*line(~par,4:6);
point(:,~par) = edge.p1(:,~par) + repmat(t(~par),3,1) .* line.directorVector(:,~par);

% set points outside of edge to [NaN NaN NaN]
point(:,t<0) = NaN;
point(:,t>1) = NaN;

indices = point(1,:)==point(1,:);

end
