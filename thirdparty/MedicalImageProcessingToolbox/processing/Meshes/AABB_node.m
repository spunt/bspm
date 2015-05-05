classdef AABB_node < handle
    % This class defines a AABB_node
    % by Alberto Gomez, 2012
    % inspired by:
    %   http://www.kreationsedge.com/?page_id=27
    %
    %  André Borrmann, Stefanie Schraufstetter, Christoph van Treeck, Ernst Rank, "An iterative, octree-based algorithm for distance computation
    % between polyhedra with complex surfaces"
    
    properties(GetAccess = 'public', SetAccess = 'public')
        % cell with the son nodes
        leftNodes=[];
        rightNodes=[];
        depth;
        mesh = MeshType();
        bounds = [0 0 0 0 0 0]';
        index_among_broders;
        
    end
    
    properties(GetAccess = 'public', SetAccess = 'private')
        unique_triangles_id=[];
        
    end
    
    methods(Access = public)
        %constructor
        function obj = AABB_node(mesh,depth,unique_triangles_id,index_among_broders)
            if (nargin > 0)
                obj.mesh= mesh;
                obj.bounds([1 3 5]) = min(mesh.points);
                obj.bounds([2 4 6]) = max(mesh.points);
                obj.depth = depth;
                obj.unique_triangles_id = unique_triangles_id;
                obj.index_among_broders = index_among_broders;
            end
        end
        
        
        
        function index = FindBestAxis(obj)
            %divide this box into two boxes - pick a better axis
            iAxisResult = [0 0 0]; % stores how close end result is, the lower the better
            %center = (obj.bounds([1 3 5])+obj.bounds([2 4 6]))/2;
            pointsOfTriangles = cat(3,obj.mesh.points(obj.mesh.triangles(:,1),:),...
                obj.mesh.points(obj.mesh.triangles(:,2),:),...
                obj.mesh.points(obj.mesh.triangles(:,3),:));
            faceCenter = mean(pointsOfTriangles,3);
            center =  mean(faceCenter);
            
            
            for iAxis=1:3
                
                left = numel(nonzeros(faceCenter(:,iAxis) <= center(iAxis)));
                right = size(obj.mesh.triangles,1)- left;
                iAxisResult(iAxis) = abs(left-right);
            end
            
            [dummy index] = min(iAxisResult);
        end
        
        
        function theyIntersect = intersect_ray_with_triangle(obj, points, normal)
            % Says whether a ray intersects a triangle. The node hasto have
            % just one triangle
            % Triangle to ray implemented following:
            % T. Moller and B. Trumbore, "Fast, Minimum Storage
            % Ray/Triangle Intersection", Journal of Graphics Tools,
            % 2(1):21--28, 1997.
            %compute BB of ray
            npts = size(points,1);
            theyIntersect = ones(npts,1);
            test_cull=false;
            
            if obj.mesh.ntriangles > 1
                disp('WARNING: this should never be reached');
                return;
            end
            
            epsilon =1E-06;
            c_points = obj.mesh.points(obj.mesh.triangles(1,:),:); % points of the triangle
            % c_normal = cross(c_points(2,:)-c_points(1,:),c_points(3,:)-c_points(1,:));
            
            edge1 = c_points(2,:)-c_points(1,:);
            edge2 = c_points(3,:)-c_points(1,:);
            
            pvec = cross(normal,edge2);
            dt = dot(edge1,pvec);
            
            if test_cull
                %------------
                if  dt < epsilon
                    %theyIntersect=theyIntersect*0;% Ray is parallel to the triangle so it will never intersect
                    theyIntersect=[];% Ray is parallel to the triangle so it will never intersect
                    return;
                end
                
                tvec = points-ones(npts,1)*c_points(1,:);
                u= dot(tvec, ones(npts,1)*pvec,2);
                
                theyIntersect(u<0 | u>dt)=0;
                
                qvec = cross(tvec,ones(npts,1)*edge1,2);
                v =dot ( ones(npts,1)*normal, qvec,2);
                
                theyIntersect(v<0 | (v+u)>dt)=0;
                
                
                 t = (edge2 * qvec');
                 inv_det = 1/dt;
                 t = t*inv_det;
                 u = u*inv_det;
                 v = v*inv_det;
                
                
                %-----------
            else
                
                if dt>-epsilon && dt < epsilon
                    %theyIntersect=theyIntersect*0;% Ray is parallel to the triangle so it will never intersect
                    theyIntersect=[];
                    return;
                end
                
                inv_det = 1/dt;
                
                tvec = points-ones(npts,1)*c_points(1,:);
                u = dot(tvec, ones(npts,1)*pvec,2)*inv_det;
                
                theyIntersect(u<=0 | u>=1)=0;
                
                qvec = cross(tvec,ones(npts,1)*edge1,2);
                v =dot ( ones(npts,1)*normal, qvec,2)*inv_det;
                
                theyIntersect(v<=0 | (v+u)>=1)=0;
                
                
                 t = (edge2 * qvec')*inv_det;
                 theyIntersect (t<=0)=0; % Select triangles only in the direction of the line.
            end
            
            theyIntersect = find(theyIntersect);
            
        end
        
        function theyIntersect = intersect(obj,otherNode)
            % Check if AABBs fail to overlap in any direction
            theyIntersect=true;
            
            noOverlap = obj.bounds([1 3 5])>otherNode.bounds([2 4 6]) | ...
                obj.bounds([2 4 6])<otherNode.bounds([1 3 5]);
            
            if sum(noOverlap)
                theyIntersect=false;
            end
            
        end
        
        function theyIntersect = intersectBB(obj,bb)
            % Check if AABBs fail to overlap in any direction
            
            
            n=size(bb,2); % number of points
            
            
            theyIntersect = ~(bb(2,:) < obj.bounds(1) | bb(4,:) < obj.bounds(3) | bb(6,:) < obj.bounds(5) | ...
                bb(1,:) > obj.bounds(2) | bb(3,:) > obj.bounds(4) | bb(5,:) > obj.bounds(6) );
            
            theyIntersect = theyIntersect(:);
            
            
            
            
        end
        
        function intersection = computeIntersection(obj,otherNode)
            % This function returns the intersection between two nodes. A
            % this stage, it is required that each node has only one
            % triangle. The intersection will be meaningfull only if the
            % two triangles intersect at a segment, oterwise NaN will be
            % returned
            %
            %   The <= comparisons make that we ignore points, which is a
            %   feature and not a bug
            %
            % TBD: orientation of the segment may help to distinguish
            % between union, intersection and difference. Think of this.
            
            intersection = [NaN NaN NaN; NaN NaN NaN];
            
            %triangles vertices
            c_points = obj.mesh.points(obj.mesh.triangles(1,:),:);
            o_points = otherNode.mesh.points(otherNode.mesh.triangles(1,:),:);
            
            % 1. Test if one triangle is totally contained at one side of
            % the plane defined by the other. If this is the case, return.
            
            c_normal = cross(c_points(2,:)-c_points(1,:),c_points(3,:)-c_points(1,:));
            c_normal = c_normal/norm(c_normal);
            
            o_signed_distances =  c_normal * (o_points - c_points([1 1 1],:))'; % distances from points in o to c
            if (nnz(o_signed_distances<=0)==0 || nnz(o_signed_distances<=0)==3)
                return;
            end
            
            % Now test the other triangle
            o_normal = cross(o_points(2,:)-o_points(1,:),o_points(3,:)-o_points(1,:));
            o_normal = o_normal/norm(o_normal);
            
            c_signed_distances =  o_normal * (c_points - o_points([1 1 1],:))'; % distances from points in c to o
            if (nnz(c_signed_distances<=0)==0 || nnz(c_signed_distances<=0)==3)
                return;
            end
            % The triangles may overlap in the line which intersects
            % planes.
            [line_p line_n] = intersectionPlanePlane(c_normal',c_points(1,:)',o_normal',o_points(1,:)');
            if (isnan(line_p))
                return;
            end
            % project the vertices of the first triangle onto the line
            c_projectedPoints = (c_points - [1 1 1]'*line_p)*line_n;
            o_projectedPoints = (o_points - [1 1 1]'*line_p)*line_n;
            %compute the parameter t
            
            if (mean(sign(c_signed_distances))<0)
                c_two = c_signed_distances<0;
                c_one = c_signed_distances>0;
            else
                c_two = c_signed_distances>0;
                c_one = c_signed_distances<0;
            end
            
            if (mean(sign(o_signed_distances))<0)
                o_two = o_signed_distances<0;
                o_one = o_signed_distances>0;
            else
                o_two = o_signed_distances>0;
                o_one = o_signed_distances<0;
            end
            
            c_t =  sort(c_projectedPoints(c_one)*[1 ; 1]  - (c_projectedPoints(c_one)*[1 ; 1]-c_projectedPoints(c_two)).*(c_signed_distances(c_one) * [1 ; 1])./(c_signed_distances(c_one)*[1 ; 1]-c_signed_distances(c_two)'));
            o_t =  sort(o_projectedPoints(o_one)*[1 ; 1]  - (o_projectedPoints(o_one)*[1 ; 1]-o_projectedPoints(o_two)).*(o_signed_distances(o_one) * [1 ; 1])./(o_signed_distances(o_one)*[1 ; 1]-o_signed_distances(o_two)'));
            
            % now find out whether there is an intersection and of which
            % type
            
            [ leftmost_t leftmost_triangle ] = max([ c_t(1) o_t(1)]);
            [ rightmost_t rightmost_triangle ] = min([ c_t(2) o_t(2)]);
            
            if (leftmost_t > rightmost_t)
                return;
            end
            if (leftmost_t == rightmost_t)
                intersection = line_p + leftmost_t*line_n';
                % intersection is just one point!!
                return;
            end
            
            % get the intersection points
            intersection = [1 ; 1]*line_p + [leftmost_t ; rightmost_t]*line_n';
            
        end
        
    end
    
    
    
    
end

