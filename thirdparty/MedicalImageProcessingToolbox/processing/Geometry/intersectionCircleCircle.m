function [point, theta1, theta2] = intersectionCircleCircle(circle1,circle2)
% function point = intersectionCircleCircle(circle1, circle2)
% This function calculates the intersection point (if any) between two
% circles in 3D
%

% inhtersection between circle planes:

%% SOME ERROR HERE: THE INTERSECTION LINE-CIRCLE SEEMS TO BE CORRECTLY 
%% CALCULATED BUT THE ANGLES SEEMS NOT. CHECK GRAPHICALLY FOR THE THETA 
%% ANGLE OF THE CIRCLE RANGE AND FOR THE CALCULATED THETA

d_th = 1e-06; % maximum distance accepted as same point
angle_th = 1e-03; % 1 degree tolerance

[pout, nout] = intersectionPlanePlane(circle1.direction(:,3),circle1.centre,...
    circle2.direction(:,3),circle2.centre);

point = NaN;
theta1 = NaN;
theta2 = NaN;

if nnz(isnan(nout)) || nnz(isnan(pout))
    return
end

if 0
    
    m1 = circle3DMesh(circle1.centre,circle1.radius,circle1.direction,'resolution',35,'theta',circle1.theta);
    m1b = circle3DMesh(circle1.centre,circle1.radius,circle1.direction,'resolution',35,'theta',[-180 180]);
    m2 = circle3DMesh(circle2.centre,circle2.radius,circle2.direction,'resolution',35,'theta',circle2.theta);
    m2b = circle3DMesh(circle2.centre,circle2.radius,circle2.direction,'resolution',35,'theta',[-180 180]);
    figure; axis equal;
    plotpoints(m1b.points','r--');
    hold on;  plotpoints(m2b.points','b--'); hold off;
    hold on;  plotpoints(m1.points','r-','LineWidth',3); hold off;
    hold on;  plotpoints(m2.points','b-','LineWidth',3); hold off;
    hold on; plotPlane(circle1.centre,circle1.direction(:,3),'scale',30,'color',[1 0 0],'opacity',0.5); hold off
    hold on; plotPlane(circle2.centre,circle2.direction(:,3),'scale',30,'color',[0 0 1],'opacity',0.5); hold off
    hold on; plotpoints([pout-nout*140 pout+nout*100]','k-','LineWidth',2); hold off;
    hold on; plotAxis( circle1.direction, circle1.centre,'scale',10); hold off
    hold on; plotAxis( circle2.direction, circle2.centre,'scale',10); hold off
     hold on;  plotpoints(intersection_point_wc1(1:3,:)','g.','MarkerSize',30); hold off;
     
end


%% find the points where the line intersects the circles - point 1
M1 = [circle1.direction circle1.centre; 0 0 0 1];
line_point_circleCoords1 = M1 \ [pout pout+nout; 1 1];

v = line_point_circleCoords1(1:2,2)-line_point_circleCoords1(1:2,1);
d = norm(v); 
D = det(line_point_circleCoords1([1 2],:));

disc = circle1.radius^2*d^2-D^2;

if disc<0
    return
end

intersection_point1(1,:) = [ (D*v(2)+ sign(v(2))*v(1)*sqrt(disc))/d^2   (D*v(2)- sign(v(2))*v(1)*sqrt(disc))/d^2];
intersection_point1(2,:) = [ (-D*v(1)+ abs(v(2))*sqrt(disc))/d^2   (-D*v(1)- abs(v(2))*sqrt(disc))/d^2 ];
intersection_point_wc1 = M1 * [intersection_point1 ; 0 0 ; 1 1];

v = v/d;

if 0
   % display 2D world 
   % something wrong here
   n=100;
   angl = (0:n-1)/(n-1)*2*pi;
   x = circle1.radius*cos(angl);
   y = circle1.radius*sin(angl);
   angl2 = (circle1.theta(1):(circle1.theta(2)-circle1.theta(1))/(n-1):circle1.theta(2))*pi/180;
   x2 = circle1.radius*cos(angl2);
   y2 = circle1.radius*sin(angl2);
   figure, 
   plot(x,y,'b--');
   hold on;
   plot(x2,y2,'b-','LineWidth',3);
   plot(circle1.radius*[0 cos(circle1.theta(1)*pi/180)],circle1.radius*[0 sin(circle1.theta(1)*pi/180)],'r-');
   plot(circle1.radius*[0 cos(circle1.theta(2)*pi/180)],circle1.radius*[0 sin(circle1.theta(2)*pi/180)],'r-');
   % line
   plotpoints([line_point_circleCoords1(1:3,1)+100*[v;0] line_point_circleCoords1(1:3,1)-10*[v;0]]','.-k','MarkerSize',20);
   plot(circle1.radius*[0 cos(theta_1(2)*pi/180)],circle1.radius*[0 sin(theta_1(2)*pi/180)],'b--','LineWidth',2);
   %point
   theta_1 = mod(atan2(intersection_point1(2,:),intersection_point1(1,:)),2*pi)*180/pi;
   plot(circle1.radius*cos(theta_1(2)*pi/180),circle1.radius*sin(theta_1(2)*pi/180),'r.','MarkerSize',15);
   plot(intersection_point1(1,:),intersection_point1(2,:),'go','MarkerSize',5); % KEY THING!
   hold off;
   axis equal
   xlabel('x')
   
end
%% THE GRAPH ABOVE SHOWS THAT THE INTERSECTION POINTS ARE BADLY COMPUTED!!

%hold on; plotpoints(intersection_point_wc1(1:3,:)','.m','MarkerSize',20); hold off

%% find the points where the line intersects the circles - point 2
M2 = [circle2.direction circle2.centre; 0 0 0 1];
line_point_circleCoords2 = M2 \ [pout pout+nout; 1 1];
v2 = line_point_circleCoords2(1:2,2)-line_point_circleCoords2(1:2,1);
v2 = v2([1 2]);
d2 = norm(v2);
D2 = det(line_point_circleCoords2([1 2],:));

disc2 = circle2.radius^2*d2^2-D2^2;

if disc2<0
    return
end

intersection_point2 = [(D2*v2(2)+ sign(v2(2))*v2(1)*sqrt(disc2))/d^2 (D2*v2(2)- sign(v2(2))*v2(1)*sqrt(disc2))/d2^2
    (-D2*v2(1)+ abs(v2(2))*sqrt(disc2))/d2^2    (-D2*v2(1)- abs(v2(2))*sqrt(disc2))/d2^2];
intersection_point_wc2 = M2 * [intersection_point2 ; 0 0 ; 1 1];

v2 = v2/d;
%hold on; plotpoints(intersection_point_wc2(1:3,1)','.m','MarkerSize',20); hold off

% aligned points (1 with 1 and 2 with 2)
distances = sum((intersection_point_wc2-intersection_point_wc1).^2);
samepoint = find(distances<d_th);
point_1_order = [1 2];
if ~numel(samepoint)
    % try with anti-aligned points: 1 with 2 and 2 with 1
    point_1_order = [2 1];
    distances = sum((intersection_point_wc2-intersection_point_wc1(:,point_1_order)).^2);
    samepoint = find(distances<d_th);
    if ~numel(samepoint)
        return;
    end
end

%% see if the point in circle coordinates is within the range
theta_1 = mod(atan2(intersection_point1(2,:),intersection_point1(1,:))*180/pi,360);
theta_2 = mod(atan2(intersection_point2(2,:),intersection_point2(1,:))*180/pi,360);

point_index1  = angleInRange( theta_1,circle1.theta, angle_th);
point_index2  = angleInRange( theta_2,circle2.theta, angle_th);


if numel(point_index1)==0 || numel(point_index2 )==0 || ~numel(intersect(point_index2,point_1_order(point_index1) )) || ...
        (numel(samepoint)==1 && point_index1~=samepoint)
    return;
end

point_index = intersect(point_index2,point_1_order(point_index1) );

% if numel(point_index)>1
%     a=1;
% end


point = intersection_point_wc2(1:3,point_index);
theta1 = theta_1(point_index1);
theta2 = theta_2(point_index);

end