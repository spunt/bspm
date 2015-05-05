classdef QuadTetMeshType < handle 
    % This class defines a trinagulated mesh
    % by Alberto Gomez, 2011
    %
    
    properties(GetAccess = 'public', SetAccess = 'public')
        npoints=0;
        ncells=0;
        points=[];% npointsx3 matrix
        cells=[];%ntrianglesx3 matrix
        bounds = [0 0 0 0 0 0];
        attributes = AttributeType(); % for labels, luts, etc
    end
    
    methods(Access = public)
        %constructor
        function obj = QuadMeshType(npoints,ncells)
            if(nargin > 0)
                if(nargin ==1 )
                    % copy constructor
                    obj.npoints = npoints.npoints;
                    obj.ncells = npoints.ncells;
                    obj.points = npoints.points;
                    obj.cells = npoints.cells;
                    for i=1:numel(npoints.attributes)
                        obj.attributes(i)=AttributeType(npoints.attributes(i));
                    end
                    
                else
                    obj.npoints = npoints;
                    obj.ncells = ncells;
                    obj.points = zeros(npoints,3);
                    obj.cells = zeros(ncells,10); % this is quadratic meshes!
                    obj.attributes(1)=AttributeType(obj.npoints);
                end
            end
        end
        
                
     
  
        
        
        function p = find_attribute(obj,at_name)
            % returns the index of the attribute if found, -1 otherwise.
            n_attributes = numel(obj.attributes);
            p=-1;
            for i=1:n_attributes
                if (strcmp(at_name, obj.attributes(i).name))
                    p = i;
                    return;
                end
            end
            
            
        end
    
        % Cleanup the points
        function removeDuplicatedPoints(obj)
            % This function removes duplicated points from the mesh and updates topology accordingly
            [nodup_points, m, n] = unique(obj.points,'rows');
            clear nodup_points;
            m_sorted  = sort(m);
            obj.points = obj.points(m_sorted ,:);
            old_m_sorted = 1:obj.npoints;
            points_to_be_removed = setdiff(old_m_sorted , m_sorted');
            replacements = m(n(points_to_be_removed))';
            
            % Remove duplicated labels
            for i=1:numel(obj.attributes)
                if (size(obj.attributes(i).attribute_array,1)==obj.npoints)
                    obj.attributes(i).attribute_array =obj.attributes(i).attribute_array(m_sorted,:);
                end
            end
            
            
            % Replace indices with duplicateds
            for i=1:numel(points_to_be_removed)
                obj.cells(obj.cells == points_to_be_removed(i))=replacements(i);
            end
            
            % remove points and update topology
            
            while(numel(points_to_be_removed))
                obj.cells(obj.cells>points_to_be_removed(1))=obj.cells(obj.cells>points_to_be_removed(1))-1;
                points_to_be_removed(1)=[];
                points_to_be_removed=points_to_be_removed-1;
            end
            
            obj.npoints = size(obj.points,1);
            
        end
        
     
        
     
        
    end % methods
    
    
    
    
end

function n = cyclicNext(a,b)
n = a+1;
if (n>b)
    n = 1;
end
end

