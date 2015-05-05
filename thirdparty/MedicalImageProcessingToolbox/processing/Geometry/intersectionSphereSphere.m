function out = intersectionSphereSphere(sphere1,sphere2)
% function point = intersectionSphereSphere(sphere1,sphere2)
%
% Computes the intersection between the two spheres. 
% the inputs are two sphere structures of the following shape:
% sphere.radius = r;
% sphere.centre = [x y z]';
%
% The output is 
%   out.direction = eye(3);
%   out.centre = [x y z]';
%   out.radius = r;
%   out.sphere(1) = sphere1
%   out.sphere(2) = sphere2
        
a = norm(sphere2.centre-sphere1.centre);

x = (a^2+sphere1.radius^2-sphere2.radius^2)/(2*a);



normal = (sphere2.centre-sphere1.centre);
normal = normal/norm(normal);

doIntersect = abs(x)<sphere1.radius;
if ~doIntersect
    out.direction = NaN(3,3);
    out.centre = [NaN NaN NaN]';
    out.radius = NaN;
    return;
end

P = sphere1.centre+x*normal;
r = sqrt(sphere1.radius^2-x^2);
 
out.direction =eye(3);
[xv, yv ]= vtkMathPerpendiculars(normal,pi/2);
out.direction(:,1) =xv;
out.direction(:,2) =yv;
out.direction(:,3) =normal;
out.centre = P;
out.radius = r;
%out.sphere(1)=sphere1;
%out.sphere(2)=sphere2;

end