function out = intersectionLineSphere(line,sphere)
% function point = intersectionLineSphere(line,sphere2)
%
% line is:
%   line.normal
%   line.point
% sphere is
%   sphere.radius
%   sphere.centre
%
% out is the two points or NaN if there is no intersection
        

tmp =  (line.normal' * (line.point-sphere.centre))^2 -...
    (line.point-sphere.centre)'*(line.point-sphere.centre) +...
    sphere.radius^2;

if tmp<0
    out = NaN(3,2);
    return
end

d1 = line.normal' * (line.point-sphere.centre) +  ...
    sqrt(tmp );

d2 = line.normal' * (line.point-sphere.centre) -  ...
    sqrt( tmp  );


point1 = line.point - d1*line.normal;
point2 = line.point - d2*line.normal;


out = [point1 point2];

end