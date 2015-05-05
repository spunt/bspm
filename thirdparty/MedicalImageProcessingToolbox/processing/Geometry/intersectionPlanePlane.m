function [pout, nout] = intersectionPlanePlane(n1,p1,n2,p2)
% returns a point and the director vector which define the line of
% intersection of two planes
%
% Inputs:
%   n1,n2 normal of planes;
%   p1,p2 points of planes
%
%   

pout = NaN;
nout=NaN;

nout = cross(n1,n2);
if norm(nout)<1e-10
    return;
end
nout = nout/norm(nout);



M = [n1  n2]';
b = [ n1'*p1 ; n2'*p2];

pout = M\b;

% % plot
% 
% m1  = planeMesh(p1, n1,'scale',1000);
% m2  = planeMesh(p2, n2,'scale',1000);
% 
% figure; 
% hold on;
% viewMesh(m1,'color',[0 0 1]);
% viewMesh(m2,'color',[1 0 0]);
% hold on; plotpoints(pout','g*','MarkerSize',60);hold off;
% hold off;


% if ( (n1(1)~=0 && n2(1)~=0) || (n1(2)~=0 && n2(2)~=0))
%     M = [ n1(1) n1(2) ; n2(1) n2(2)];
%     if (abs(det(M))<10^-14)
%         return;
%     end
%     
%     p = M\b;
%     pout(3)=0;
%     pout(1:2)=p;
%     
% elseif ( (n1(1)~=0 && n2(1)~=0) || (n1(3)~=0 && n2(3)~=0))
%     
%     M = [ n1(1) n1(3) ; n2(1) n2(3)];
%     if (abs(det(M))<10^-14)
%         return;
%     end
%     
%     p = M\b;
%     pout(2)=0;
%     pout([1 3])=p;
%     
% else
%     M = [ n1(2) n1(3) ; n2(2) n2(3)];
%     if (abs(det(M))<10^-14)
%         return;
%     end
%     
%     p = M\b;
%     pout(1)=0;
%     pout([2 3])=p;
%     
% end



end