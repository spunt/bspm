function [theta, point] = intersectionCircleTriangle(circle,triangle)
% function point = intersectionCircleTriangle(circle,triangle)
% triangle is provided as three points (column vectors)
% must be assumed that triangle and circle are in the same plane

Mcirc = [ circle.direction circle.centre; 0 0 0 1];
triangle_2D = Mcirc \ [triangle; 1 1 1];


% obtain frustum planes
%% get the 4 corners of the base of the frustum
side= sqrt(frustum.mesh.npoints/2);
fourcorners = frustum.mesh.points([1 side  side*side (side-1)*side+1 ]+frustum.mesh.npoints/2,:);
ath = -.5; % 2 degrees of tolerance
points = fourcorners - ones(4,1)*frustum.centre';
% project along z direction
vectors_magnitude = sqrt(sum(points.^2,2));
vectors = points./[vectors_magnitude  vectors_magnitude  vectors_magnitude ];
%plotpoints(vectors,'*')

%% calculate the 4 normal vectors

frustum_face_vector(1,:) = cross(vectors(1,:),vectors(2,:));
frustum_face_vector(2,:) = cross(vectors(2,:),vectors(3,:));
frustum_face_vector(3,:) = cross(vectors(3,:),vectors(4,:));
frustum_face_vector(4,:) = cross(vectors(4,:),vectors(1,:));


% sort the vectors by the direction
frustum_face_vector_projected = frustum.directions\frustum_face_vector';
angles_with_x = mod(atan2(frustum_face_vector_projected(2,:),frustum_face_vector_projected(1,:))*180/pi,360);

%Always take as first the one over 0 but closest

[~, faces_sorted] = sort( angles_with_x);

% the planes are now defined by 'points' and frustum_face_vector (row vectors and points)
% the intersection between a plane and a circle is

% 1. obtain the intersection point in w c.


point = [];
excess_point=[];
whichface= [];
thetafrustum = [];
phifrustum = [];

for face=1:4
    % 1.1 Intersection between frustum plane and circle plane (a line)
    
    
    line.normal = cross(frustum_face_vector(face,:), circle.direction(:,3))';
    
    m = [frustum_face_vector(face,:)'  circle.direction(:,3)]';
    b = [frustum_face_vector(face,:)*fourcorners(face,:)' ; circle.direction(:,3)'*circle.centre(:)];
    line.point = m\b;
    
    % project both the line and the circle into the circle plane to do
    % the calculation in 2D
    
    line2D.normal = [line.normal'*circle.direction(:,1) line.normal'*circle.direction(:,2)];
    line2D.normal = line2D.normal /norm(line2D.normal );
    M = [circle.direction circle.centre; 0 0 0 1];
    line2D.point  = M\[ line.point(:) ; 1];
    line2D.point = line2D.point(1:2);
    
    circle2D.centre = M\[ circle.centre(:) ; 1];
    circle2D.centre = circle2D.centre(1:2); % this HAS TO be (0,0)
    circle2D.radius = circle.radius;
    
    % intersect 2D line with 2D circle
    
    
    D = det([line2D.point(:) line2D.point(:)+line2D.normal(:)]);
    laplacian = circle2D.radius^2-D^2;
    if laplacian <=0
        % no intersection or just tangent
        intersection{face}=NaN(3,2);
        continue;
    end
    
    intersection2D(1,1)=D*line2D.normal(2)+sign(line2D.normal(2))*line2D.normal(1)*sqrt(laplacian); % x comp point 1
    intersection2D(1,2)=D*line2D.normal(2)-sign(line2D.normal(2))*line2D.normal(1)*sqrt(laplacian);% x comp point 2
    intersection2D(2,1)=-D*line2D.normal(1)+abs(line2D.normal(2))*sqrt(laplacian);% y comp point 1
    intersection2D(2,2)=-D*line2D.normal(1)-abs(line2D.normal(2))*sqrt(laplacian);% y comp point 1
    
    tmp = M*[intersection2D; 0 0; 1 1];
    intersection{face}=tmp(1:3,:);
    
    
    % Now to check if the point is inside. For that, check the
    % angles in both directions, thats the only way!!
    
    points_in_frustum_coordinates = [frustum.directions frustum.centre(:) ; 0 0 0 1]\[ intersection{face}; 1 1];
    points_in_frustum_coordinates=points_in_frustum_coordinates(1:3,:);
    % Check if the point is inside the frustum. First compute radius
    norm_points_in_frustum_coordinates=sqrt(sum(points_in_frustum_coordinates.^2 ));
    in_range_radius = find(~sum(sign(frustum.radiuses(:)*[1 1]-[1 1]'*norm_points_in_frustum_coordinates)));
    %  if ~numel(in_range_radius )
    %      continue;
    %  end
    % Check if the point is inside the frustum. Second compute angles
    vectors = points_in_frustum_coordinates./...
        [norm_points_in_frustum_coordinates;
        norm_points_in_frustum_coordinates;
        norm_points_in_frustum_coordinates];
    
    % TODO: SAVE ALL POINTS
    % The get only the two which look better within range
    %     angles_limit = frustum.angles+[-ath ath -ath ath]/2;
    %
    
    if 0
        
        mesh2 = MeshType(frustum.mesh) ;
        mesh2.points =  [frustum.directions frustum.centre(:) ; 0 0 0 1]\ [ frustum.mesh.points ones(frustum.mesh.npoints,1)]';
        mesh2.points = mesh2.points(1:3,:)';
        figure, viewMesh(mesh2,'opacity',0.5);
        hold on; plotpoints(points_in_frustum_coordinates(:,2)','g.','MarkerSize',15); hold off
        axis on;
        grid on;
        xlabel('x')
        ylabel('y')
        hold on; plotAxis(eye(3),[0 0 0]','scale',25); hold off;
    end
    
    
    angles_limit = frustum.angles;
    angles = [atan2(vectors(1,:),vectors(3,:))
        atan2(vectors(2,:),vectors(3,:))]*180/pi; % first row: angle1, second row: angle2
    
    excess(1,:)= angles(1,:)-angles_limit(1)*[1 1]; % excess_thetaMin
    excess(2,:) = angles_limit(2)*[1 1]-angles(1,:); % excess_thetaMax
    excess(3,:) = angles(2,:)-angles_limit(3)*[1 1]; % excess_phiMin
    excess(4,:) = angles_limit(4)*[1 1]-angles(2,:); % excess_phiMax
    
    excessNeg = excess;
    excessNeg(excessNeg >0)=0;
    excess = sum(excessNeg);
    
    excess_point = [excess_point excess];
    point = [point intersection{face}];
    whichface= [whichface faces_sorted(face)*ones(1,2)];
    thetafrustum = [thetafrustum  angles(1,:)];
    phifrustum = [phifrustum angles(2,:)];
    
end



if ~numel(point)
    point = NaN(3,2);
    theta=[NaN NaN];
    return;
end
% 2. translate point into theta coordinates relative to the
M = [ circle.direction circle.centre; 0 0 0 1];
points_in_circle_coordinates = M\[ point; ones(1,size(point,2))];

% circle axis
% calculate angles

theta=atan2(points_in_circle_coordinates(2,:),points_in_circle_coordinates(1,:))*180/pi;


% select the two points which better fit
[~,candidate_points] = sort(excess_point);
% check for duplicates (This might happen when the points are very close to a corner)
there_is_duplicate = true;
while there_is_duplicate
    dist_between_points = norm(point(:,candidate_points(end))-point(:,candidate_points(end-1)));
    if dist_between_points/circle.radius < 0.05
        % points are the same
        excess_point(candidate_points(end-1))=-Inf;
        [~,candidate_points] = sort(excess_point);
    else
        there_is_duplicate = false;
    end
    
end
candidate_points = candidate_points(end-1:end);

% check that these two points are inside the frustum
if excess_point(candidate_points(1)) < ath
    point = NaN(3,2);
    theta=[NaN NaN];
    return;
end

point = point(:,candidate_points);
whichface = whichface(candidate_points);
thetafrustum = [thetafrustum(candidate_points);phifrustum(candidate_points)];
theta = theta(candidate_points);

 % get the point between the two theta values and check if it is inside the
 % frustum
[~,midrange]=angleInRange(0,theta,1E-03);
midrange=[midrange midrange+180];
mid_point_in_circle_coordinates = circle.radius*[cos(midrange*pi/180) ;sin(midrange*pi/180) ;0 0];
mid_point = M*[mid_point_in_circle_coordinates;1 1];
mid_point_in_frustum_coordinates = [frustum.directions frustum.centre(:) ; 0 0 0 1]\mid_point;

mid_point_frustum_angles= [atan2(mid_point_in_frustum_coordinates(1,:),mid_point_in_frustum_coordinates(3,:))
                           atan2(mid_point_in_frustum_coordinates(2,:),mid_point_in_frustum_coordinates(3,:))]*180/pi;
a1 = angleInRange(mid_point_frustum_angles(1,:),frustum.angles(1:2),1E-03);
a2 = angleInRange(mid_point_frustum_angles(2,:),frustum.angles(3:4),1E-03);

point_inside = intersect(a1,a2);

i = [1 2 1];
i = i(point_inside:point_inside+1);

theta = theta(i);
point = point(:,i);
thetafrustum = thetafrustum(:,i); % first row: angle1. Second row: angle 2
whichface=whichface(i);

end