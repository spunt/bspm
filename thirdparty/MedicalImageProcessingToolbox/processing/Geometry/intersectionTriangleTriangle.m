function [points triangleIndex] = intersectionTriangleTriangle(t1,t2)
% returns the two points of intersection between two triangular patches and
% the triangles they belong two (1 or 2)
%
% t1 and t2 are a 3 x 3 matrix where each column is a point counter clockwise, of the
% triangle
%


% TODO modify this to process many triangles at the same time ... ?
n1 = cross(t1(:,1)-t1(:,2),t1(:,1)-t1(:,3));
n1 = n1/norm(n1);

n2 = cross(t2(:,1)-t2(:,2),t2(:,1)-t2(:,3));
n2 = n2/norm(n2);



[pout nout] = intersectionPlanePlane(n1,t1(:,1),n2,t2(:,1));

% find the two points where this line intersects with the triangle t1

normals_line(:,1)= (t1(:,2)-t1(:,1))'/norm(t1(:,2)-t1(:,1));
normals_line(:,2)= (t1(:,3)-t1(:,2))'/norm(t1(:,3)-t1(:,2));
normals_line(:,3)= (t1(:,1)-t1(:,3))'/norm(t1(:,1)-t1(:,3));

normals_line(:,4)= (t2(:,2)-t2(:,1))'/norm(t2(:,2)-t2(:,1));
normals_line(:,5)= (t2(:,3)-t2(:,2))'/norm(t2(:,3)-t2(:,2));
normals_line(:,6)= (t2(:,1)-t2(:,3))'/norm(t2(:,1)-t2(:,3));

t = [t1 t2];

points=[];
triangleIndex=[];
for i=1:6
    % intersection between two lines
    Up =  normals_line(:,i);
    Uq =  nout(:);
    
    P0 = t(:,i);
    Q0 = pout(:);
    
    current_point =intersectionLineLine(Up,P0,Uq,Q0); % the closest point between both
    
    if (current_point(1)~=current_point(1))
        continue;
    end
    
    % see if the point is in the opposite triangle
    
    % if the point is inside both triangles, add it
    % first change coordinates
    if (i<4) % trinagle 1
        [x y] = vtkMathPerpendiculars(n2,pi/2);
        M = [x y n2];
        t_2D = M\ t2; % triangle always has 3 points
    else
        [x y] = vtkMathPerpendiculars(n1,pi/2);
        M = [x y n1];
        t_2D = M\ t1; % triangle always has 3 points
    end
    
    current_point_2D = M\current_point;
    
    if ( PointInTriangle(current_point_2D([1 2]),t_2D([1 2],1),t_2D([1 2],2),t_2D([1 2],3)))
        points = [points current_point];
        triangleIndex(size(points,2)) = (i>3)+1;
    end
    
    
end

    
end




