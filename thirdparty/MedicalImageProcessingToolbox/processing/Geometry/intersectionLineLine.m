function point = intersectionLineLine(Up,P0,Uq,Q0)
% function point = intersectionLineLine(n1,p1,n2,p2)
%
% Computes the closest point between two lines
% Must be column points
        
        b = -1*[(P0-Q0)'*Up; (P0-Q0)'*Uq];
        A = [ Up'*Up -Uq'*Up;  Up'*Uq -Uq'*Uq];
        if (abs(det(A))<10^(-10))
            point = [NaN NaN NaN]';
        else
            lambda = A\b;

            P1 = P0 + lambda(1)*Up;
            Q1 = Q0 + lambda(2)*Uq;
            point = (P1+Q1)/2;
        end
        
end