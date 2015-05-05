function [out,midrange,rangeWidth] = angleInRange( angle,range, angle_th)
% This function gives true if the angle(s) (in degrees) is in the given
% range. Angles must be given in mod 360

while range(2)<range(1)
    range(2) =range(2)+360;
end
rangeWidth = range(2)-range(1);
midrange = mean(range);
out = zeros(size(angle));
for an =numel(angle):-1:1
    angle(an) = angle(an)-360;
    while angle(an)<(range(1)-angle_th)
        angle(an)=angle(an)+360;
    end
    out(an)=angle(an)<=(range(2)+angle_th);    
end


out=find(out);


end

