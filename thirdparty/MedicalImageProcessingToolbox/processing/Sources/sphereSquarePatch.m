function out_mesh  = sphereSquarePatch(c,r,varargin)
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


resolution = 16; % half resolution for phi than for theta
theta = [-180 180].*pi/180;
phi = [ -180 180].*pi/180;
sphere_axes = eye(3);
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
       elseif  (strcmp(varargin{i},'direction'))
           sphere_axes= varargin{i+1}; 
           i = i+2;
       end
   end
end

         
    
    anglex_span = (theta(2)-theta(1))/resolution;
    angley_span = (phi(2)-phi(1))/resolution;
            
    angles_x = theta(1):anglex_span:theta(2);
    angles_y = phi(1):angley_span:phi(2);
    sp_counter = 0;
    points = zeros(numel(angles_x)*numel(angles_y),3);
    for ax = angles_x
        for ay = angles_y
            sp_counter = sp_counter+1;
          %  points(sp_counter,:) = r*[sin(ax) cos(ax)*sin(ay) cos(ax)*cos(ay) ];
          %points(sp_counter,:) = r*[sin(ax)*cos(ay) cos(ax)*sin(ay) cos(ax)*cos(ay) ];
          points(sp_counter,:) = r*[sin(ax)*cos(ay) cos(ax)*sin(ay) sqrt(1-cos(ax).^2*sin(ay).^2-cos(ay).^2*sin(ax).^2) ];
        end
    end
            
    % topology of the surfaces
            
    points_surfaces_2D = points(:,1:2);
    DT = DelaunayTri(points_surfaces_2D);
    triangles = DT.Triangulation;
            
    %  indices of points in the sides of the surfacein  counter-clockwise order from the beam source
    indices_of_surfaces_sides_up = [ 1:numel(angles_y):((numel(angles_x)-1)*numel(angles_y)+1) ...
                ((numel(angles_x)-1)*numel(angles_y)+1):1:(numel(angles_x)*numel(angles_y)) ...
                (numel(angles_x)*numel(angles_y)):-numel(angles_y):numel(angles_x) ...
                numel(angles_x):-1:1];
            indices_of_surfaces_sides_up(numel(angles_y))=[];
            indices_of_surfaces_sides_up(numel(angles_y)+numel(angles_x)-1)=[];
            indices_of_surfaces_sides_up(2*numel(angles_y)+numel(angles_x)-2)=[];
            
            
            out_mesh = MeshType(size(points,1),size(triangles,1));
            out_mesh.triangles = triangles;
            % the points have to be transformed to the appropiate orientation and origin
            
            
            
            
          out_mesh.points=( [ sphere_axes c(:); 0 0 0 1] *[points' ; ones(1,size(points,1))])';
            
            

   
   
  

end
