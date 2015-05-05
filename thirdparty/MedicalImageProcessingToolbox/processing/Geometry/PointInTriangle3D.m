function flag = PointInTriangle3D(p, a,b,c)
% returns true if the point p is inside the triangle defined by a,b,c (in 3D)
% p, a b and c are column points
    

    flag = SameSide(p,a, b,c) .* SameSide(p,b, a,c) .* SameSide(p,c, a,b);
        
end

function flag = SameSide(p1,p2, a,b)
% This function and the next one are taken from http://www.blackpawn.com/texts/pointinpoly/default.html 

    npts = size(p1,2);
    p1 = [p1; zeros(1,npts)];
    %p2Large = p2*ones(1,npts);
    
    aLarge = [a; 0]*ones(1,npts);
    bLarge = [b; 0]*ones(1,npts);
    
    cp1 = cross(bLarge-aLarge, p1-aLarge,1);
    cp2 = cross([b-a; 0], [p2; 0]-[a; 0],1);
    cp2Large = cp2*ones(1,npts);
    
    epsilon = 10^-10;
    
    flag = dot(cp1,cp2Large,1) >= (0-epsilon);
        
end


