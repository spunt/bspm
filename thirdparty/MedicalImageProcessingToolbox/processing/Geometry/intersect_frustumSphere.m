function [intersection_points theta_circle theta_frustum circles]= intersect_frustumSphere(frustum, sphere_radius)
% function point = intersectionCircleFrustum(n1,p1,n2,p2)
%

% obtain frustum planes
%% get the 4 corners of the base of the frustum
side= sqrt(frustum.mesh.npoints/2);
fourcorners = frustum.mesh.points([1 side  side*side (side-1)*side+1 ]+frustum.mesh.npoints/2,:);
points = fourcorners - ones(4,1)*frustum.centre';
% project along z direction

vectors_magnitude = sqrt(sum(points.^2,2));
vectors = points./[vectors_magnitude  vectors_magnitude  vectors_magnitude ];
%plotpoints(vectors,'*')

%% calculate the 4 normal vectors

frustum_face_vector(1,:) = cross(vectors(1,:),vectors(2,:)); % vectors are sorted counter-clockwise (seen from the frustum centre)
frustum_face_vector(2,:) = cross(vectors(2,:),vectors(3,:));
frustum_face_vector(3,:) = cross(vectors(3,:),vectors(4,:));
frustum_face_vector(4,:) = cross(vectors(4,:),vectors(1,:));


% sort the vectors by the direction
frustum_face_vector_projected = frustum.directions\frustum_face_vector';
angles_with_x = atan2(frustum_face_vector_projected(2,:),frustum_face_vector_projected(1,:))*180/pi;
angles_with_x(angles_with_x<0)=angles_with_x(angles_with_x<0)+360;
%Always take as first the one over 0 but closest

[~, faces_sorted] = sort( angles_with_x);

%frustum_face_vector = frustum_face_vector(faces_sorted,:);

% the planes are now defined by 'points' and frustum_face_vector (row vectors and points)
% the intersection between a plane and a circle is

% 1. obtain the intersection point in w c.
intersection_points = zeros(3,4);
theta_circle = zeros(2,4); % theta coordinate for each point (one per circle that produces that point)
theta_frustum= zeros(2,4);

face_index = [1 2 3 4 1 ];
for face=1:4
    % intersection between two contiguous frustum faces => line
    % [~, nout] =...
    %     intersectionPlanePlane(frustum_face_vector(face_index(face),:),fourcorners(face_index(face),:),...
    %     frustum_face_vector(face_index(face+1),:),fourcorners(face_index(face+1),:));
    
    
    
    intersection_points(:,face) = frustum.centre + vectors(face_index(face),:)'*sphere_radius;
    
    circles(face).direction =eye(3);
    [xv yv ]= vtkMathPerpendiculars(frustum_face_vector(face,:)',pi/2);
    circles(face).direction(:,1) =xv;
    circles(face).direction(:,2) =yv;
    circles(face).direction(:,3) =frustum_face_vector(face,:);
    
    circles(face).centre =   frustum.centre;
    circles(face).radius = sphere_radius;
    
end
%% now calculate angles (theta)

fsorted = faces_sorted([4 1 2 3]);

for face = 1:4
    % calculate angles in circle coordinates
    
    M = [ circles(face).direction circles(face).centre; 0 0 0 1];
    points_in_circle_coordinates = M\[ intersection_points(:,[face_index(face) face_index(face+1)]); ones(1,2)];
    
    theta=mod((atan2(points_in_circle_coordinates(2,:),points_in_circle_coordinates(1,:))*180/pi),360); % first row is from the first side, second row from the second side
    
      if 0
         figure; axis equal
         m1 = circle3DMesh(circles(face).centre,circles(face).radius,circles(face).direction,'resolution',35,'theta',[-180 180]);
         m1b = circle3DMesh(circles(face).centre,circles(face).radius,circles(face).direction,'resolution',35,'theta',theta);
         hold on;  plotpoints(m1.points','r--','LineWidth',1); hold off;
         hold on;  plotpoints(m1b.points','r-','LineWidth',3); hold off;
         hold on; plotAxis( circles(face).direction, circles(face).centre,'scale',10); hold off
         hold on; plotpoints(intersection_points(:,[face_index(face) ])','b.','MarkerSize',25); hold off
     end
    
    
    theta_circle(1,face_index(face))=theta(1);
    theta_circle(2,face_index(face+1))=theta(2);

    circles(face).theta = theta;
    circles(face).points = fsorted([face_index(face) face_index(face+1)]);
    %mesh  = circle3DMesh(circles(face).centre,circles(face).radius,circles(face).direction,'resolution',35,'theta',circles(face).theta);
    %figure, hold on; plotpoints(mesh.points','k-','LineWidth',3); hold off; axis equal
    %mesh2  = circle3DMesh(circles(face).centre,circles(face).radius,circles(face).direction,'resolution',35,'theta',[-150 150]);
    % hold on; plotpoints(mesh2.points','k--','LineWidth',1); hold off; 
    %hold on; plotAxis( circles(face).direction*Mrot', circles(face).centre,'scale',10); hold off
    % calculate angles in frustum coordinates
    M = [frustum.directions frustum.centre(:) ; 0 0 0 1];
    points_in_frustum_coordinates = M\[ intersection_points(:,face_index(face)); 1];
    points_in_frustum_coordinates=points_in_frustum_coordinates(1:3,:);
    norm_points_in_frustum_coordinates=sqrt(sum(points_in_frustum_coordinates.^2 ));
    vecs = points_in_frustum_coordinates./[norm_points_in_frustum_coordinates;norm_points_in_frustum_coordinates;norm_points_in_frustum_coordinates];
    theta_frustum(:,face) = [atan2(vecs(1,:),vecs(3,:));atan2(vecs(2,:),vecs(3,:))]*180/pi; % first row: angle1, second row: angle2
    
%      if 0
%         
%         m1 = circle3DMesh(circles(face).centre,circles(face).radius,circles(face).direction,'resolution',35,'theta',circles(face).theta);
%         hold on;  plotpoints(m1.points','r-','LineWidth',3); hold off;
%     end
     
    
end

intersection_points = intersection_points(:,fsorted);
circles = circles(faces_sorted);
theta_circle = theta_circle(:,faces_sorted);
theta_frustum = theta_frustum(:,faces_sorted);



end