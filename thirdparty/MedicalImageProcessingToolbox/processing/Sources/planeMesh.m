function m  = planeMesh(point, normalVector,varargin)
% m =  planeMesh(point, normalVector)
% m =  planeMesh(point, normalVector,options)
%
% creates a MeshType object of a sphere of radius r and center c
%
% options       meaning     default
% -------       -------     -------
%
% 's'           scale       1


s = 1; % half resolution for phi than for theta


% Argument reading

if (size(varargin,2)>0)
    i = 1;
    while (i <= size(varargin,2))
        if(strcmp( varargin{i} , 'scale'))
                s = varargin{i+1};
        end
        i = i+1;
    end
end

   
   % Create plane
   
   [x y ]   =vtkMathPerpendiculars(normalVector,pi/2);
   
  
   points(:,1)=point(:)-s/2*x-s/2*y;
   points(:,2)=point(:)-s/2*x+s/2*y;
   points(:,4)=point(:)+s/2*x-s/2*y;
   points(:,3)=point(:)+s/2*x+s/2*y;
      
   % Create conectivity
   
        
   polygons = [1 2 3; 1 3 4];
   
   
     
   % create mesh
   m = MeshType(size(points,1),size(polygons,1)); 
   m.points =points';
   m.triangles = polygons;
   
   
  
end

