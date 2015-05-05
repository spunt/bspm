function m  = sphereMesh(c,r,varargin)
% m = sphereMesh(c,r)
% m = sphereMesh(c,r,options)
%
% creates a MeshType object of a sphere of radius r and center c
%
% options       meaning     default
% -------       -------     -------
%
% 'resolution'      sphereRes   16
% 'theta'           angle       [0 360]
% 'phi'             angle       [0 180]
% 'rotatez'         angle       0
% 'rotatey'         angle       0
% 'rotatex'         angle       0
% 'rotationOrder'   int         {0} = Rx - Ry - Rz
%                               1 = Ry - Rx - Rz

global closed;
resolution = 8; % half resolution for phi than for theta
theta = [0 360].*pi/180;
phi = [0 180].*pi/180;
rotatez = 0;
rotatey = 0;
rotatex = 0;
rotationOrder=0;
nOptions = size(varargin,2);
% Argument reading

if (nOptions > 0)
    if (rem(nOptions,2)~=0)
        disp('Error in the arguments, please check the option list');
        return;
    end
    i=1;
   while(i<=nOptions)
       
       if (strcmp(varargin{i},'resolution'))
           resolution = varargin{i+1};
           i = i+2;
       elseif  (strcmp(varargin{i},'theta'))
           theta = varargin{i+1}*pi/180;
           i = i+2;
       elseif  (strcmp(varargin{i},'phi'))
           phi = varargin{i+1}*pi/180; 
           i = i+2;
       elseif  (strcmp(varargin{i},'rotatez'))
           rotatez = varargin{i+1}*pi/180; 
           i = i+2;
       elseif  (strcmp(varargin{i},'rotatey'))
           rotatey = varargin{i+1}*pi/180; 
           i = i+2;
       elseif  (strcmp(varargin{i},'rotatex'))
           rotatex = varargin{i+1}*pi/180; 
           i = i+2;
        elseif  (strcmp(varargin{i},'rotationOrder'))
           rotationOrder= varargin{i+1}; 
           i = i+2;
       end
   end
end

phi(1)=max(phi(1),0);
phi(2)=min(phi(2),180);
theta(1)=max(theta(1),0);
theta(2)=min(theta(2),360);

closed = true;
if (theta(2)-2*pi< theta(1))
    closed = false;
end
    
   % --------------------------
   
   thetaResolution = resolution;
   phiResolution = resolution;
   
   
   numPts = phiResolution * resolution + 2; % phires * thetares +2
   numPolys = phiResolution *2* resolution;
   
   
   
   
   % Create sphere
   %t = (-thetaResolution:2:thetaResolution)/thetaResolution*pi;
   deltaTheta = (theta(2)-theta(1))/thetaResolution;
   deltaPhi = (phi(2)-phi(1))/phiResolution;
   t = theta(1):deltaTheta:theta(2);
   %p = (-phiResolution:2:phiResolution)'/phiResolution*pi/2;
   p = (phi(1):deltaPhi:phi(2))';

      
   % Create conectivity
   
   polygons = [];
   points = [];
   
   numpoles = 0;
   poleNorthFound = false;
   poleSouthFound = false;
  
   for phi_=p';
       for theta_=t;       
           % north pole
           if (~poleNorthFound && ~phi_)
               poleNorthFound = true;
               numpoles = numpoles+1;
               for i=1:thetaResolution
                   vertex1 = i+1;
                   vertex2 = cyclicNext((vertex1),[2 thetaResolution+numpoles]);
                   vertex3 = 1; % the pole
                   polygons = [polygons ; vertex1 vertex2 vertex3];
               end
                x = r*sin(phi_)*cos(theta_);
                y = r*sin(phi_)*sin(theta_);
                z = r*cos(phi_);
                points = [points; [x y z]+c'];
           end
           % north pole
           if (~poleSouthFound && phi_>=pi)
               poleSouthFound = true;
               numpoles = numpoles+1;
               
               for i=1:thetaResolution
                   vertex1 = size(points,1)-thetaResolution+i;
                   vertex2 = cyclicNext((vertex1),[size(points,1)+1-thetaResolution size(points,1)]);
                   vertex3 = size(points,1)+1; % the pole
                   polygons = [polygons ; vertex1 vertex2 vertex3];
               end
                x = r*sin(phi_)*cos(theta_);
                y = r*sin(phi_)*sin(theta_);
                z = r*cos(phi_);
                  points = [points; [x y z]+c'];
           end
            if (phi_ && phi_<pi && theta_ < t(end))
               x = r*sin(phi_)*cos(theta_);
                y = r*sin(phi_)*sin(theta_);
                z = r*cos(phi_);
                 points = [points; [x y z]+c'];
                           
           end
            
           
      end
   end
   
   % band conectivity
   
   for i=2:(phiResolution-1);
     for j=1:(thetaResolution)
         
           % conectivity
           vertex1 = (i-2)*thetaResolution+j+1;
           vertex2 = cyclicNext(vertex1,[(i-2)*thetaResolution+2 (i-2)*thetaResolution+thetaResolution+1]);
           vertex3 = cyclicNext(vertex1+thetaResolution-1,[(i-1)*thetaResolution+2 (i-1)*thetaResolution+thetaResolution+1]);
           polygons = [polygons ; vertex1 vertex2 vertex3];
           
           vertex1 = (i-2)*thetaResolution+j+1;
           
           vertex3 = cyclicNext(vertex1+thetaResolution-1,[(i-1)*thetaResolution+2 (i-1)*thetaResolution+thetaResolution+1]);
           vertex2 = cyclicPrevious(vertex3,[(i-1)*thetaResolution+2 (i-1)*thetaResolution+thetaResolution+1]);
           polygons = [polygons ; vertex1 vertex2 vertex3];
           
           
           
      end
   end
   
   Rx = eye(3);
   Ry = eye(3);
   Rz = eye(3);
   
   if (rotatex~=0)
       Rx = [1 0 0
             0 cos(rotatex) -sin(rotatex)
             0 sin(rotatex) cos(rotatex)];
   end
   
   if (rotatey~=0)
       Ry = [cos(rotatey) 0 sin(rotatey)
           0 1 0
           -sin(rotatey) 0 cos(rotatey)];
   end
   
   if (rotatez~=0)
       Rz = [cos(rotatez) -sin(rotatez) 0
           sin(rotatez) cos(rotatez) 0
           0 0 1];
   end

   switch (rotationOrder)
       case 0
           M = Rz * Ry * Rx;
       case 1
           M = Rz * Rx * Ry;
   end
   
   for i=1:size(points,1)
           points(i,:) = (M * (points(i,:)'-c(:)) +c(:))';
   end
   

   
   
     
   % create mesh
   m = MeshType(size(points,1),size(polygons,1)); 
   m.points =points;
   m.triangles = polygons;
   
   
   % just for testing

end


   function n = cyclicNext(a,b)
   global closed;
   % a index
   % b=[b1 b2] min and max index
  % if (a<b(1))
  %     a=b(1);
  % end
  
    
    n = a+1;
    if (n>b(2))
        if (closed)
            n = b(1);
        else
            n=a;
        end
    end
    
   end

   function n = cyclicPrevious(a,b)
   global closed;
   % a index
    
    n = a-1;
    if (n<b(1))
        if (closed)
            n = b(2);
        else
            n=a;
        end
    end
    
        
    
end