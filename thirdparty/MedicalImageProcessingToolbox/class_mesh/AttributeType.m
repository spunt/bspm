classdef AttributeType
    % This class defines an attribute in a polydata file
    % by Alberto Gomez, 2011
    %
    
     properties(GetAccess = 'public', SetAccess = 'public')
         attribute;
         name;
         numComp=1; % between 1 and 4
         type;
         lookup_table;
         nelements; %number of tuples
         
         attribute_array=[];
    end
    
    methods(Access = public)
        %constructor
        function obj = AttributeType(npoints)
            obj.lookup_table='default';
            if (nargin==1 && strcmp(class(npoints),'AttributeType'))
                % copy constructor
                obj.attribute_array = npoints.attribute_array;
                obj.nelements = npoints.nelements;
                obj.name = npoints.name;
                obj.attribute=npoints.attribute;
                obj.type=npoints.type;                
                obj.numComp=npoints.numComp;
                obj.lookup_table=npoints.lookup_table;
            elseif (nargin>0)
                obj.attribute_array = zeros(npoints,1);
                obj.nelements = npoints;
                obj.name = 'default';
                obj.attribute='field';
                obj.type='float';                
            end
        end
        
  
 
       
    end % methods
    
    
        
    
end


