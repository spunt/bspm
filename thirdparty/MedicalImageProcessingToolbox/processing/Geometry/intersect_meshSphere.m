function [theta_frustum,circles,points]= intersect_meshSphere(frustum, sphere_radius, mesh)
%
%%
% from the mesh select elegible vertices to intersect. That happens only if
% at least two points of a triangle are at different sides of the spherical
% surface.

th_dist = 1E-06;
th_angle = 0.5; % 0.5 degrees
points = NaN;

% calculate angles in frustum coordinates
M = [frustum.directions frustum.centre(:) ; 0 0 0 1];
points_in_frustum_coordinates = M\[ mesh.points'; ones(1,mesh.npoints)];
points_in_frustum_coordinates=points_in_frustum_coordinates(1:3,:);
radiuses = sqrt(sum(points_in_frustum_coordinates.^2));

% calculate theta and phi angles

vecs = points_in_frustum_coordinates./[radiuses;radiuses;radiuses];
theta_frustum = atan2(vecs(1,:),vecs(3,:))*180/pi;
phi_frustum = atan2(vecs(2,:),vecs(3,:))*180/pi; % first row: angle1, second row: angle2

out_thetaMin = find(theta_frustum<(frustum.angles(1)-th_angle));
out_thetaMax = find(theta_frustum>(frustum.angles(2)+th_angle));
out_phiMin = find(phi_frustum<(frustum.angles(3)-th_angle));
out_phiMax = find(phi_frustum>(frustum.angles(4)+th_angle));

out_of_frustum = unique([out_thetaMin out_thetaMax out_phiMin out_phiMax]);
inside_frustum = ones(mesh.npoints,1);
inside_frustum(out_of_frustum) = 0;

% intersects = zeros(mesh.npoints,1);
% for i =1:mesh.ntriangles
%     sides =unique((radiuses(mesh.triangles(i,:)) - sphere_radius*ones(1,3))>0); % at least one point at either side of the sphere
%     insides = nnz(inside_frustum(mesh.triangles(i,:))); % at least one point inside the frustum
%     intersects(i)= numel(sides)>1 & insides >0;
%     % hold on; plotpoints(mesh.points(mesh.triangles(i,:),:),'b.','MarkerSize',31); hold off
%  end
% 
% intersecting_triangles  = find(intersects)
intersecting_triangles  = 1:mesh.ntriangles;

% hold on; plotpoints(mesh.points(mesh.triangles(intersecting_triangles,:),:),'r.','MarkerSize',31); hold off

if ~numel(intersecting_triangles)
    circles = NaN;
    theta_frustum = NaN;
    return
end


%% find the intersection between each triangle and the sphere, which
% must be a plane. note that  some points might well be outside the
% frustum.
circles = [];

points = [];
% if sphere_radius >52 && sphere_radius < 54
%     disp('stop here')
% end
for i = intersecting_triangles(:)'
    % the intersection of one triangle with one sphere is a circular
    % arch curve. This arch is defined by two points. these two points
    % are the intersection of the triangle sides which have the two
    % points at opposite sides of the sphere and the sphere
    
    
    
    l1p1 = mesh.points(mesh.triangles(i,2),:);
    l1p2 = mesh.points(mesh.triangles(i,1),:);
    
    l2p1 =  mesh.points(mesh.triangles(i,3),:);
    l2p2 = mesh.points(mesh.triangles(i,1),:);
    
    l3p1 =  mesh.points(mesh.triangles(i,3),:);
    l3p2 = mesh.points(mesh.triangles(i,2),:);
    
    
    line1.normal = l1p1(:)- l1p2(:);
    line1.normal = line1.normal/norm(line1.normal);
    line1.point = l1p2(:);
    
    line2.normal = l2p1(:)- l2p2(:);
    line2.normal = line2.normal/norm(line2.normal);
    line2.point = l2p2(:);
    
    line3.normal = l3p1(:)- l3p2(:);
    line3.normal = line3.normal/norm(line3.normal);
    line3.point = l3p2(:);
    
    sphere.radius = sphere_radius;
    sphere.centre= frustum.centre(:);
    
    points1 = intersectionLineSphere(line1,sphere);
    points2 = intersectionLineSphere(line2,sphere);
    points3 = intersectionLineSphere(line3,sphere);
    
    % Discard points which are outside the two vertices for each line
    
    %        projection on the line
    range1 = sort([ line1.normal'*l1p1(:)        line1.normal'*l1p2(:)]);
    pr1 = line1.normal'*points1;
    out_triangle_ind = zeros(1,6);
    out_triangle_ind(1:2) = (pr1<=range1(1) | pr1>=range1(2));
    %points1=points1(:,pr1>range1(1) & pr1<range1(2));
    if points1~=points1
        out_triangle_ind(1:2)=[1 1];
        points1  = NaN(3,2);
    end
    
    range2 = sort([ line2.normal'*l2p1(:)        line2.normal'*l2p2(:)]);
    pr2 = line2.normal'*points2;
    out_triangle_ind(3:4) = pr2<=range2(1) | pr2>=range2(2);
    %points2= points2(:,pr2>range2(1) & pr2<range2(2));
    if points2~=points2
        out_triangle_ind(3:4)=[1 1];
        points2  = NaN(3,2);
    end
    
    
    range3 = sort([ line3.normal'*l3p1(:)        line3.normal'*l3p2(:)]);
    pr3 = line3.normal'*points3;
    %points3= points3(:,pr3>range3(1) & pr3<range3(2));
    
    out_triangle_ind(5:6)=pr3<=range3(1) | pr3>=range3(2);
    
    if points3~=points3
        out_triangle_ind(5:6)=[1 1];
        points3  = NaN(3,2);
    end
    
    if ~nnz(out_triangle_ind==0)
        continue;
    end
    
    points_candidates = [points1 points2 points3];
    out_triangle = find(out_triangle_ind);
    in_triangle = setdiff(1:6,out_triangle );
    
    
    % now that we have the points in the triangle, calculate the angles in
    % circle coordinates
    
    
    circle(1).direction = eye(3,3);
    vz = cross(l1p2-l1p1,l2p2-l2p1);
    circle(1).direction(:,3) = vz/norm(vz);
    [vx,vy]= vtkMathPerpendiculars(circle(1).direction(:,3),pi/2);
    circle(1).direction(:,[1 2])=[vx vy];
    if points1(1)==points1(1)
        circle(1).centre =  ((points1(:,1)-frustum.centre)'*circle(1).direction(:,3))*circle(1).direction(:,3)+frustum.centre;
        circle(1).radius = norm(points1(:,2)-circle(1).centre );
    elseif points2(1)==points2(1)
        circle(1).centre =  ((points2(:,1)-frustum.centre)'*circle(1).direction(:,3))*circle(1).direction(:,3)+frustum.centre;
        circle(1).radius = norm(points2(:,2)-circle(1).centre );
    else
        circle(1).centre =  ((points3(:,1)-frustum.centre)'*circle(1).direction(:,3))*circle(1).direction(:,3)+frustum.centre;
        circle(1).radius = norm(points3(:,2)-circle(1).centre );
    end
    circle(1).theta = [-180 180]; % this does not harm and will be changed later on
    
    
    Mcirc = [ circle(1).direction circle(1).centre; 0 0 0 1];
    triangle_2D = Mcirc \ [points_candidates(:,in_triangle); ones(1,numel(in_triangle))];
    v2D = Mcirc \ [ mesh.points(mesh.triangles(i,:),:)'; 1 1 1]; % vertices of the triangle in 2D
    theta_in_triangle = atan2(triangle_2D(2,:),triangle_2D(1,:))*180/pi;
    
  
    % 1. calculate the intersection of the circle with the frustum
    %[theta_in_frustum, point_in_frustum, whichface, thetafrustum] = intersectionCircleFrustum(circle(1),frustum);
    theta_in_frustum = intersectionCircleFrustum(circle(1),frustum);
    
    % 2. calculate the intersection of the circle with the triangle
    % NOTE : there can be more than two points in the triangle (rare but
    % possible):
    % 3 points: one side of the triangle cut the sphere and another is
    % tangent
    % 4 points: Two sides of the triangle cut the sphere
    % 5 points: Two sides of the triangle cut the sphere and the third is
    % tangent
    % 6 points: The three sides of the triangle cut the sphere 
    
    if size(theta_in_triangle,2)==2
    
        % to verify the order of theta, lets take the point in the middle of
        % the range ad see if it is inside the triangle
        [~,midrange] = angleInRange(0,theta_in_triangle ,1E-03);
        p_midrange = circle(1).radius*[cos(midrange*pi/180) sin(midrange*pi/180) 0]';

        % calculate the intersection of the circle with the frustum
        arch_that_overlaps=1;
        range_thetas = 1;
        sorted_theta_in_triangle = theta_in_triangle;
    elseif size(theta_in_triangle,2)==3
        disp('3 intersections')
    elseif size(theta_in_triangle,2)==4
        [sorted_theta_in_triangle,sorting_order] = sort(theta_in_triangle);
        [~,midrange1] = angleInRange(0,sorted_theta_in_triangle(1:2) ,1E-03);
         p_midrange(:,1) = circle(1).radius*[cos(midrange1*pi/180) sin(midrange1*pi/180) 0]';
        [~,midrange2] = angleInRange(0,sorted_theta_in_triangle(3:4) ,1E-03);
         p_midrange(:,2) = circle(1).radius*[cos(midrange2*pi/180) sin(midrange2*pi/180) 0]';
        
        arch_that_overlaps=[1 1];      
        range_thetas = [1 2];
        
        % create the two circles
        
        
    elseif size(theta_in_triangle,2)==5
        disp('5 intersections')
    elseif size(theta_in_triangle,2)==6
        disp('6 intersections')
    end
        % see if there is any overlap between the arc in the frustum and the
        % arc in the triangle
    clear final_theta_in_triangle;
  
    % convert point to 3D
    for ii=range_thetas
        p_midrange_3D = Mcirc*[ p_midrange(:,ii) ; 1];
        if ~PointInTriangle(p_midrange(1:2,ii), v2D(1:2,1),v2D(1:2,2),v2D(1:2,3)) %|| ... %if point is outside the triangle
           % ~PointInFrustum(p_midrange_3D,frustum)  % or point is outside the frustum
            final_theta_in_triangle(ii,:) = sorted_theta_in_triangle(([ii+1/2 ii]-1)*2+1);
        else
            final_theta_in_triangle(ii,:) = sorted_theta_in_triangle(([ii ii+1/2]-1)*2+1);
        end     
        if ~numel(angleInRange(final_theta_in_triangle(ii,:),theta_in_frustum,1E-03)) &&...
           ~numel(angleInRange(theta_in_frustum,final_theta_in_triangle(ii,:),1E-03))
            arch_that_overlaps(ii)=0;
        end
    end

    if ~nnz(arch_that_overlaps)
        %disp('Does not really intersect')
        continue;
    end
        overlapping_indices = find(arch_that_overlaps);
        clear final_theta_sorted circle_;
        
    for ii=1:numel(overlapping_indices ) % for each of these we will produce one circle
        
        % there is some overlap. we take the most restrictive and we're done!
        if ~numel(angleInRange(final_theta_in_triangle(overlapping_indices(ii),:),theta_in_frustum,1E-03))
            % all interval theta2 is contained within
            % theta1
            sorting_order = [1 3 4 2];
        elseif ~numel(angleInRange(theta_in_frustum,final_theta_in_triangle(overlapping_indices(ii),:),1E-03))
            % all interval theta1 is contained within
            % theta2
            sorting_order = [3 1 2 4];
        else
            % Both intervals do partial overlap
            if numel(angleInRange(final_theta_in_triangle(overlapping_indices(ii),2),theta_in_frustum,1E-03)) ||...
                    numel(angleInRange(theta_in_frustum(1),final_theta_in_triangle(overlapping_indices(ii),:),1E-03))
                sorting_order = [1 3 2 4];
            elseif numel(angleInRange(final_theta_in_triangle(overlapping_indices(ii),1),theta_in_frustum,1E-03)) || ...
                    numel(angleInRange(theta_in_frustum(2),final_theta_in_triangle(overlapping_indices(ii),:),1E-03))
                sorting_order = [3 1 4 2];
            else
                disp('ERROR when ordering theta in intersectionBetweenSpheres.m')
            end
        end
        final_theta_sorted(ii,:) = [final_theta_in_triangle(overlapping_indices(ii),:) theta_in_frustum];
        final_theta_sorted(ii,:) =final_theta_sorted (ii,sorting_order);
    
        circle_(ii,1).direction = circle(1).direction;
        circle_(ii,1).centre = circle(1).centre;
        circle_(ii,1).radius = circle(1).radius;
        % produce the points in circle coordinates
        final_points_2D = [circle_(ii,1).radius*cos(final_theta_sorted(ii,[2 3])*pi/180) ; circle_(ii,1).radius*sin(final_theta_sorted(ii,[2 3])*pi/180) ; 0 0; 1 1];
        final_points = Mcirc *final_points_2D;
        
        % assign the theta
        circle_(ii,1).theta = final_theta_sorted(ii,[2 3]);

        % assign the points
        % check if any point is already in
        isAlreadyUsed = [0 0];

        if numel(points)
            for iii=[1 2]
                dist = sum((points-final_points(1:3,iii)*ones(1,size(points,2))).^2);
                previous_equal = find(dist<th_dist);
                if numel(previous_equal)
                    [~,previous_equal]=min(dist); % in case there are several below th_dist
                    
                    isAlreadyUsed(iii)=previous_equal;
                end
            end
        end
        point_index = [size(points,2)+1  size(points,2)+2-(isAlreadyUsed(1)>0)];
        points = [points final_points(1:3,~isAlreadyUsed)];

        point_index(isAlreadyUsed>0) = isAlreadyUsed(isAlreadyUsed>0) ;
        circle_(ii,1).points = point_index;
    end
    
    if 0
        % plot spherical patch
        mesh_patch = sphereSquarePatch(frustum.centre,sphere_radius,'resolution',16,'theta',frustum.angles([1 2]),'phi',frustum.angles([3 4]),'direction',frustum.directions);
        figure
        hold on;
        viewMesh(mesh_patch,'opacity',0.5,'color',[0.5 0 0]);
        hold off;
        %plot triangle
        hold on;
        plotpoints(mesh.points(mesh.triangles(i,[1 2 3 1]),:));
       % plotpoints(mesh.points(mesh.triangles(i,find(sides)),:),'k*');
        hold off
        
        hold on; plotpoints(final_points(1:3,:)','g.','MarkerSize',25); hold off
        hold on; plotCircleArch(circle_,'Color',[0 1 0]); hold off
       % hold on; viewMesh(mesh, 'triangles',intersecting_triangles); hold off;
    end
    
    
    
    circles = [circles ;circle_];
end

if ~numel(circles)
    points = NaN;
    theta_frustum=NaN;
end

end