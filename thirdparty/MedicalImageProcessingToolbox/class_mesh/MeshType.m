classdef MeshType < handle
    % This class defines a trinagulated mesh
    % by Alberto Gomez, 2011
    %
    
    properties(GetAccess = 'public', SetAccess = 'public')
        npoints=0;
        ntriangles=0;
        points=[];% npointsx3 matrix
        triangles=[];%ntrianglesx3 matrix
        bounds = [0 0 0 0 0 0];
        attributes = AttributeType(); % for labels, luts, etc
    end
    
    methods(Access = public)
        %constructor
        function obj = MeshType(npoints,ntriangles)
            if(nargin > 0)
                if(nargin ==1 )
                    % copy constructor
                    obj.npoints = npoints.npoints;
                    obj.ntriangles = npoints.ntriangles;
                    obj.points = npoints.points;
                    obj.triangles = npoints.triangles;
                    for i=1:numel(npoints.attributes)
                        obj.attributes(i)=AttributeType(npoints.attributes(i));
                    end
                    
                else
                    obj.npoints = npoints;
                    obj.ntriangles = ntriangles;
                    obj.points = zeros(npoints,3);
                    obj.triangles = zeros(ntriangles,3);
                    obj.attributes(1)=AttributeType(obj.npoints);
                end
            end
        end
        
        function b = ComputeBounds(obj)
            b([1 3 5])=min(obj.points,[],1);
            b([2 4 6])=max(obj.points,[],1);
            obj.bounds = b;
            
        end
        
        function v = GetVolume(obj)
            normals = obj.GetNormalAtFaces(1:obj.ntriangles);
            [~,areas]= obj.GetTriangleNormal(1:obj.ntriangles);
            p1 = obj.points(obj.triangles(:,1),:);
            p2 = obj.points(obj.triangles(:,2),:);
            p3 = obj.points(obj.triangles(:,3),:);
            pts = cat(3,p1,p2,p3);
            barycentres = mean(pts ,3);
            
            v = nansum(dot(barycentres,normals,2).*areas)/3;
            
        end
        
        
        function P = GetVertexCoordinates(obj, index)
            % get the position of a vertex at index "index", an integer
            
            P = obj.points(index,:);
        end
        
        function P  = GetTriangleVertices(obj, index)
            % get the vertices composing the triangle at the index
            P = obj.triangles(index,:);
        end
        
        function [P, area] = GetTriangleNormal(obj,index)
            % get the normal to the surface of a triangle
            
            triang = obj.triangles(index,:);
            v1 = obj.points(triang(:,2),:) - obj.points(triang(:,1),:);
            v2 = obj.points(triang(:,3),:) - obj.points(triang(:,1),:);
            cp = cross(v1,v2);
            area = sqrt(sum(cp.^2,2))/2;
            
            P = cp./([area area area]*2);
            %P = P(:)';
            
        end
        %%
        
        
        function out = laplacianSmoothing(obj)
            ignoreOrientation = true;
            out = MeshType(obj);
            h = waitbar(0,'Applying laplacian smooth');
            for i =1:obj.npoints
                nn = obj.GetVertexNeighbours(i);
                positions = obj.points(nn,:);
                out.points(i,:)=mean(positions);
                waitbar(i/obj.npoints);
            end
            close(h);
            
        end
        
        function normal = GetNormalAtVertex(obj,index, varargin)
            % get normal to a vertex
            % normal is computed as linear combination of the normals to
            % faces.
            %
            % the weigting factor for the average can be chosen by the
            % user, Gouraud being default:
            % options :
            % mode = 'Gouraud' (w=1)
            %        'area' (w=area)
            
            
            index_of_att =  obj.find_attribute('normalVectors');
            if index_of_att >0
                normal= obj.attributes(index_of_att).attribute_array(index,:);
                return
            end
            
            normal = [0 0 0];
            ignoreOrientation=false;
            center=[]; % a point inside (convex geometry) to be able to define a positicve orientation
            mode = 'Gouraud'; % w = 1
            
            for i=1:size(varargin,2)
                if (strcmp(varargin{i},'mode'))
                    mode =varargin{i+1};
                elseif (strcmp(varargin{i},'ignoreOrientation'))
                    ignoreOrientation = varargin{i+1};
                elseif (strcmp(varargin{i},'insidePoint'))
                    center= varargin{i+1};
                end
                
            end
            
            
            nn = obj.GetVertexNeighboursSorted(index,ignoreOrientation);
            
            w = ones(size(nn,1),1); % weightings
            % compute normals to the neighbor faces
            normals = zeros(length(nn),3);
            areas = zeros(length(nn),1);
            for i = 1:length(nn)
                v1 = obj.points(nn(i),:)-obj.points(index,:);
                v2 = obj.points(nn(cyclicNext(i,length(nn))),:)-obj.points(index,:);
                areas(i)=norm(cross(v1,v2))/2;
                normals(i,:) = cross(v1,v2)/areas(i)*2;
                
            end
            
            if (strcmp(mode,'Gouraud'))
                % w is all ones
                normal = sum((w*ones(1,3)).*normals,1)/sum(w);
                normal = normal/norm(normal);
            elseif (strcmp(mode,'area'))
                normal = sum((areas*ones(1,3)).*normals,1)/sum(areas);
                normal = normal/norm(normal);
            else
                normal = normals(1,:);
                normal = normal/norm(normal);
                %display('Using default mode: Gouraud');
            end
            
            if (numel(center))
                vector = obj.points(index,:)-center;
                if vector*normal'<0
                    normal=-1*normal;
                end
            end
            
        end
        %%
        function normal = GetNormalAtFaces(obj,index)
            % get normal to  faces.
            
            % get first edge
            v1 = obj.points(obj.triangles(index,2),:)-obj.points(obj.triangles(index,1),:);
            v2 = obj.points(obj.triangles(index,3),:)-obj.points(obj.triangles(index,1),:);
            
            normal = cross(v1,v2,2);
            
            % normalize result
            magnitude = sqrt(normal(:,1).^2 + normal(:,2).^2 + normal(:,3).^2);
            
            for i=1:3
                normal(:,i) = normal(:,i)./magnitude;
            end
            
            
        end
        %%
        function P = GetForwardVertexNeighbours(obj, index)
            % Returns neighbours as oriented subcession of 2 neighbours
            
            [r, c]=find(obj.triangles==index);
            extended_triangles = [obj.triangles(r,:) obj.triangles(r,:)];
            P=[];
            for i=1:numel(r)
                c = find(extended_triangles(i,:)==index,1,'first');
                P =[P; extended_triangles(i,(c+1):(c+2))];
            end
        end
        %%
        function P = GetBackwardVertexNeighbours(obj, index)
            % Returns neighbours as oriented subcession of 2 neighbours
            
            [r, c]=find(obj.triangles==index);
            extended_triangles = [obj.triangles(r,:) obj.triangles(r,:)];
            P=[];
            for i=1:numel(r)
                c = find(extended_triangles(i,:)==index,1,'last');
                P =[P; extended_triangles(i,(c-1):(c-2))];
            end
        end
        %%
        function P = GetVertexNeighbours(obj,index)
            % get a column vector with the neighbours to a vertex (excludes the vertex itself)
            [r, c]=find(obj.triangles==index);
            tri = obj.triangles(r,:);
            tri_reshaped = unique(reshape(tri,3*size(r,1),1)); % remove repeated elements
            
            P = zeros(length(tri_reshaped)-1,1);
            i=1;
            for j=1:length(tri_reshaped)
                if (tri_reshaped(j)~=index)
                    P(i)=tri_reshaped(j);
                    i = i+1;
                end
                
            end
            
        end
        %%
        function sorted_neighbours = GetVertexNeighboursSorted(obj,index, varargin)
            % get a column vector with the neighbours to a vertex (excludes the vertex itself)
            % the neighbours are given in couter-clockwise order
            
            
            ignoreOrientation=false;
            
            if  (size(varargin,2)>0)
                ignoreOrientation = varargin{1};
            end
            
            [r_, c_]=find(obj.triangles==index);
            triangles = obj.triangles(r_,:);
            trianglesExtended = [triangles triangles];
            
            
            
            
            all_neighbours = unique(triangles);
            all_neighbours(all_neighbours ==index)=[];
            added_neighbours = 0;
            current_neighbour = trianglesExtended(1,find(trianglesExtended(1,:)==index ,1,'first')+1); % we start courter clockwise
            sorted_neighbours = zeros(numel(all_neighbours),1);
            added_neighbours=added_neighbours+1;
            sorted_neighbours(added_neighbours)=current_neighbour;
            
            while(numel(nonzeros(sorted_neighbours))<numel(all_neighbours))
                [r_, c_] = find(triangles==current_neighbour);
                candidate_found=false;
                for i = 1:numel(r_)
                    c_2 = find(trianglesExtended(r_(i),:)==current_neighbour,1,'first');
                    if (ignoreOrientation)
                        candidate = [ trianglesExtended(r_(i),c_2+1) trianglesExtended(r_(i),c_2+2)];
                        candidate(candidate==index)=[];
                    else
                        candidate = trianglesExtended(r_(i),c_2+1);
                    end
                    if numel(candidate) && ~numel(find(sorted_neighbours==candidate)) && candidate~=index
                        current_neighbour=candidate;
                        added_neighbours=added_neighbours+1;
                        sorted_neighbours(added_neighbours)=current_neighbour;
                        candidate_found=true;
                        break;
                    end
                    
                end
                if ~candidate_found
                    disp('No candidate found. Check the surface orientation')
                    sorted_neighbours=[];
                    return
                end
            end
        end
        %%
        function ring = Get1Ring(obj,i)
            % i is the vertex index for which the ring is required
            
            % If there is only one triangle, this is a corner point, return
            % its triangle
            [r_, c_]=find(obj.triangles==i);
            if numel(r_)==1
                ring =r_;
                return;
            end
            
            
            % Now we assume the triangle is ordered
            ring=r_(1);
            r_(1)=[];
            
            while(numel(r_))
                vertices = obj.triangles(ring(end),[1 2 3 1 2]);
                centre_idx=find(vertices==i,1,'first');
                % Add next triangle
                next_vertex=vertices(centre_idx+2);
                
                
                triangles = obj.triangles(r_,:);
                next_triangle = find(sum(triangles==i | triangles==next_vertex,2)==2,1,'first');
                if numel(next_triangle)
                    ring = [ring r_(next_triangle)];
                    r_(next_triangle)=[];
                end
                
                % Add previous triangle
                vertices = obj.triangles(ring(1),[1 2 3 1 2]);
                centre_idx=find(vertices==i,1,'first');
                prev_vertex=vertices(centre_idx+1);
                
                triangles = obj.triangles(r_,:);
                prev_triangle = find(sum(triangles==i | triangles==prev_vertex,2)==2,1,'first');
                if numel(prev_triangle)
                    ring = [ r_(prev_triangle) ring];
                    r_(prev_triangle)=[];
                end
                
                
            end
            
        end
        
        %%
        function e = edges(obj)
            % This function returns a listof the edges
            tmp = [obj.triangles(:,[1 2]) ; obj.triangles(:,[2 3]) ; obj.triangles(:,[3 1])];
            % remove duplicates
            e=unique(sort(tmp,2),'rows');
        end
        
        %%
        function [ring,tri_sorted] = GetRingOrdered(obj,i)
            % This function returns the one ring required by the paper:
            %  M. Desbrun, M. Meyer, and P. Alliez, "Intrinsic Parametrization of
            %   Surface Meshes", Eurographics 2002
            % Arguments:
            %   i - central vertex
            % The ring is ordered so that every triangle starts with the i
            % vertex
            % This function will also correct for the orientation
            
            [r_, c_]=find(obj.triangles==i);
            
            ring=r_(1);
            r_(1)=[];
            vertices = obj.triangles(ring(end),[1 2 3 1 2]);
            tri_sorted = vertices([c_(1) c_(1)+1 c_(1)+2]);
            c_(1)=[];
            
            while(numel(r_))
                vertices = obj.triangles(ring(end),[1 2 3 1 2]);
                centre_idx=find(vertices==i,1,'first');
                % Add next triangle
                next_vertex=vertices(centre_idx+2);
                
                added=0;
                
                triangles = obj.triangles(r_,:);
                next_triangle = find(sum(triangles==i | triangles==next_vertex,2)==2,1,'first');
                if numel(next_triangle)
                    ring = [ring r_(next_triangle)];
                    vertices_ = obj.triangles(r_(next_triangle),[1 2 3 1 2]);
                    r_(next_triangle)=[];
                    tri_sorted  = [tri_sorted ; vertices_([c_(next_triangle)  c_(next_triangle)+1 c_(next_triangle)+2]) ];
                    c_(next_triangle)=[];
                    added=added+1;
                end
                
                % Add previous triangle
                vertices = obj.triangles(ring(1),[1 2 3 1 2]);
                centre_idx=find(vertices==i,1,'first');
                prev_vertex=vertices(centre_idx+1);
                
                triangles = obj.triangles(r_,:);
                prev_triangle = find(sum(triangles==i | triangles==prev_vertex,2)==2,1,'first');
                if numel(prev_triangle)
                    ring = [ r_(prev_triangle) ring];
                    vertices_ = obj.triangles(r_(prev_triangle),[1 2 3 1 2]);
                    r_(prev_triangle)=[];
                    tri_sorted  = [ vertices_([c_(prev_triangle) c_(prev_triangle)+1 c_(prev_triangle)+2]) ; tri_sorted];
                    c_(prev_triangle)=[];
                    added=added+1;
                end
                
                if added==0
                    % one of the remaining triangles is badly ordered
                    disp(['WARNING: The input mesh does not have consistent orientation at node ' num2str(i) ])
                elseif added==2
                    % Check that orientation is good
                    %disp(['WARNING: The input mesh does not have consistent orientation at node ' num2str(i) ])
                    if tri_sorted(1,3)~=tri_sorted(2,2)
                        tri_sorted(1,2)=tri_sorted(1,3);
                        tri_sorted(1,3)=tri_sorted(2,2);
                        obj.triangles(ring(1),:) = obj.triangles(ring(1),[1 3 2]);
                    end
                    
                    if tri_sorted(end,2)~=tri_sorted(end-1,3)
                        tri_sorted(end,3)=tri_sorted(end,2);
                        tri_sorted(end,2)=tri_sorted(end-1,3);
                        obj.triangles(ring(end),:) = obj.triangles(ring(end),[1 3 2]);
                    end
                end
            end
        end
        
        %%
        function out = GetValue(obj,attribute_name,positions)
            % positions must be an array of row vectors
            ind = obj.find_attribute(attribute_name);
            if (ind <= 0)
                disp([ 'Attribute ' attribute_name ' was not found '])
                out=[];
                return;
            end
            out = zeros(size(positions,1),1);
            bt = obj.ComputeBounds();
            extent = bt(2:2:end)-bt(1:2:end);
            nodim = find(extent==0);
            availabledims = setdiff(1:numel(extent),nodim);
            DT= delaunayTriangulation(obj.points(:,availabledims));
            triangle_idx = DT.pointLocation(positions);
            valid_points = find(triangle_idx==triangle_idx);
            invalues  = obj.attributes(ind).attribute_array;
            
            % out(valid_points) = invalues(valid_points); % this works fine!
            
            bary_coords = DT.cartesianToBarycentric(triangle_idx(valid_points),positions(valid_points,:));
            
            %vertices = obj.triangles(triangle_idx(valid_points),:);
            vertices = DT.ConnectivityList(triangle_idx(valid_points),:);
            
            values_at_vertices = invalues(vertices);
            
            interpolated_values = sum(values_at_vertices.*bary_coords,2);
            
            out(valid_points)=interpolated_values;
            
            
            
            
            
            
            
            
        end
        
        
        %%
        function normals=ComputeNormals(obj,varargin)
            % add normals as another attribute
            
            ignoreOrientation=false;
            center=[];
            if  (size(varargin,2)>0)
                ignoreOrientation= varargin{1};
            end
            if  (size(varargin,2)>1)
                center= varargin{2};
            end
            
            
            
            n_attributes = numel(obj.attributes);
            at = AttributeType();
            at.attribute='normals';
            at.name='normalVectors';
            at.numComp=3; % between 1 and 4
            at.type='float';
            at.nelements=obj.npoints; %number of tuples
            
            for i=1:at.nelements
                at.attribute_array(i,:) = obj.GetNormalAtVertex(i,'ignoreOrientation', ignoreOrientation,'insidePoint',center);
            end
            ind = obj.find_attribute(at.name);
            if (ind > 0)
                obj.attributes(ind)=at;
            else
                n_attributes = numel(obj.attributes);
                obj.attributes(n_attributes+1)=at;
            end
            normals = at.attribute_array;
            
        end
        
        function ComputeNormalsToFaces(obj)
            % add normals as another attribute
            
            ind = obj.find_attribute('normalVectorsFaces');
            if (ind > 0)
                return;
            end
            
            at = AttributeType();
            at.attribute='normals';
            at.name='normalVectorsFaces';
            at.numComp=3; % between 1 and 4
            at.type='float';
            at.nelements=obj.ntriangles; %number of tuples
            
            for i=1:at.nelements
                at.attribute_array(i,:) = obj.GetNormalAtFaces(i);
            end
            ind = obj.find_attribute(at.name);
            if (ind > 0)
                obj.attributes(ind)=at;
            else
                n_attributes = numel(obj.attributes);
                obj.attributes(n_attributes+1)=at;
            end
            
        end
        
        function Mesh = ExtractSurface(obj, n_attribute, value)
            % This function returns a new mesh which corresponds to the
            % sub-surface of vertex whos attribute #n_attribute have the
            % value #value
            % This function leaves the original points and attribute list
            % unchanged
            
            % 1. Check that there is the same number of points in the
            % scalar field and in the mesh
            scalar_array = obj.attributes(n_attribute).attribute_array;
            if (size(scalar_array,1)~=size(obj.points,1))
                disp('Error: associated field has too few points')
                Mesh = [];
                return;
            end
            
            index = find(scalar_array == value);
            
            % remove all triangles which contain a triangle not in index
            
            triangles_to_add=[];
            for i=1:numel(index)
                [x y]=find(obj.triangles==index(i));
                triangles_to_add = [triangles_to_add; x];
                
            end
            triangles_to_add= unique(triangles_to_add);
            tri= obj.triangles(triangles_to_add,:);
            
            Mesh = MeshType();
            Mesh.npoints = obj.npoints;
            Mesh.points = obj.points;
            Mesh.triangles = tri;
            Mesh.ntriangles = size(tri,1);
            Mesh.attributes= obj.attributes;
        end
        
        function closedMesh= closeMesh(obj)
            % this function closes the mesh by adding 1 point.
            % IMPORTANT: it will only work if there is 1 single hole
            
            closedMesh = MeshType();
            closedMesh.npoints = obj.npoints;
            closedMesh.points = obj.points;
            closedMesh.ntriangles = obj.ntriangles;
            closedMesh.triangles = obj.triangles;
            closedMesh.attributes = obj.attributes;
            
            
            [borders,borderPaths]= obj.detectBorders();
            
            nholes = numel(borderPaths);
            for i=1:nholes
                b=zeros(size(borders));
                b(unique(borderPaths{i}(:))')=1;
                bindex = find(b==1);
                borderPoints = obj.points(b==1,:);
                nborderPoints = size(borderPoints,1);
                if (nborderPoints>0)
                    %add the centroid as a new point
                    centroid = sum(borderPoints)/nborderPoints;
                    closedMesh.npoints = closedMesh.npoints+1;
                    closedMesh.points(closedMesh.npoints,:)=centroid;
                    % Build new triangles
                    
                    remainingPoints = bindex;
                    currentPoint = remainingPoints(1);
                    while (numel(remainingPoints)>1)
                        %try to look for the closest neighbor clockwise
                        minDistance = 999999999999999999; % very big number
                        for candidate=remainingPoints(:)'
                            if (currentPoint~=candidate)
                                vector = closedMesh.points(currentPoint,:)-closedMesh.points(candidate,:);
                                distance = norm(vector);
                                if (distance<minDistance )
                                    minDistance = distance;
                                    neighbor = candidate;
                                end
                            end
                        end
                        
                        %remove the currentPoint from the list
                        indexOfCurrentPoint = find(remainingPoints==currentPoint);
                        if (indexOfCurrentPoint ==1)
                            remainingPoints = remainingPoints(2:end);
                        elseif (indexOfCurrentPoint ==numel(remainingPoints))
                            remainingPoints = remainingPoints(1:(end-1));
                        else
                            remainingPoints = remainingPoints([1:(indexOfCurrentPoint-1) (indexOfCurrentPoint+1):end]);
                        end
                        
                        %build triangle
                        tri = [currentPoint neighbor closedMesh.npoints] ;
                        closedMesh.triangles = [closedMesh.triangles ; tri];
                        currentPoint = neighbor;
                    end
                    tri = [remainingPoints(1) bindex(1) closedMesh.npoints];
                    closedMesh.triangles = [closedMesh.triangles ; tri];
                    closedMesh.ntriangles = size(closedMesh.triangles,1);
                    
                end
            end
            
        end
        
        function flag = addVertex(obj,pointIn)
            % Adds a vertex (or vertices) to current mesh, updating the
            % triangulations accordingly. If there is any problem returns
            % an error code.
            % TODO: test mesh
            
            flag=1;
            % 1. Find the closest triangle
            indicesClosestTriangle = obj.findClosestTriangle(pointIn);
            
            trianglePointIndices = [ obj.triangles(indicesClosestTriangle,:) ((obj.npoints+1):(obj.npoints+size(pointIn,1)))' ];
            
            newTopology = [1 2 4 2 3 4 3 1 4]; % this will have to be reshaped
            newTriangles = trianglePointIndices(:,newTopology);
            
            % Now reshape 'newTriangles' so that from each row we get 3
            % rows.
            newTriangles = [newTriangles(:,1:3);newTriangles(:,4:6);newTriangles(:,7:9)];
            
            % Remove old triangles and add new ones
            obj.triangles(indicesClosestTriangle,:)=[];
            obj.triangles = [obj.triangles; newTriangles];
            obj.ntriangles = size(obj.triangles,1);
            
            % add new points
            obj.points = [ obj.points; pointIn];
            obj.npoints = size(obj.points,1);
            
        end
        
        function outMesh = removeVertex(obj,index,varargin)
            % Removes a vertex (or vertices) to current mesh, updating the
            % triangulations accordingly.
            %
            % This code (still to be written) can be designed to do two
            % things:
            %
            %   a) Remove the point and all triangles where this point was
            %   (this will of course produce a hole) [default]
            %
            %   b) Remove a point and re-triangulate neighbors. This will
            %   avoid holes (option 'noholes')
            %
            
            mode = 'holes'; % w = 1
            if  (size(varargin,2)>0)
                mode = varargin{1};
            end
            
            outMesh = MeshType(obj);
            outMesh.triangles = obj.triangles;
            outMesh.points = obj.points;
            
            
            if (strcmp('holes',mode))
                
                pointsToRemove = index;
                
                % remove all triangles containing any of the points to be
                % removed
                nPointsRemoved = 0;
                trianglesRemoved=[];
                while(numel(index)>0)
                    
                    [r, c]= find(outMesh.triangles==index(1));
                    outMesh.triangles(r,:)=[];
                    trianglesRemoved = [trianglesRemoved; r];
                    
                    % point per point, remove the index to the triangle list
                    
                    outMesh.points(index(1),:)=[];
                    outMesh.triangles(outMesh.triangles>index(1))=outMesh.triangles(outMesh.triangles>index(1))-1;
                    nPointsRemoved = nPointsRemoved+1;
                    index = index-1;
                    index(1)=[];
                    
                end
                
                
                % update attributes
                for i=1:numel(obj.attributes)
                    if (numel(obj.attributes(i).attribute_array)==outMesh.npoints)
                        outMesh.attributes(i).attribute_array(pointsToRemove)=[];
                    elseif (numel(obj.attributes(i).attribute_array)==outMesh.ntriangles)
                        outMesh.attributes(i).attribute_array(trianglesRemoved)=[];
                    end
                end
                
                outMesh.npoints = size(outMesh.points,1);
                outMesh.ntriangles = size(outMesh.triangles,1);
                
                
            elseif (strcmp('noholes',mode))
                
                disp(['Mode ' mode ' is not implemented yet']) ;
            else
                
                disp(['ERROR: ' mode ' is not a valid option']);
            end
            
        end
        
        
        function outMesh = removeTriangle(obj,index,varargin)
            % Removes a triangle (or triangles) to current mesh, updating the
            % triangulations accordingly.
            %
            % This code (still to be written) can be designed to do two
            % things:
            %
            %   a) Remove the triangle and leave a hole   [default]
            %
            %   b) Remove a triangle and re-triangulate neighbors. This will
            %   avoid holes (option 'noholes')
            %
            
            mode = 'holes'; % w = 1
            if  (size(varargin,2)>0)
                mode = varargin{1};
            end
            
            outMesh = MeshType(obj);
            outMesh.triangles = obj.triangles;
            outMesh.points = obj.points;
            
            if (strcmp('holes',mode))
                
                
                outMesh.triangles(index,:)=[];
                pointsToRemove = setdiff(unique(obj.triangles(:)),unique(outMesh.triangles(:)));
                
                nPointsRemoved = 0;
                while(numel(pointsToRemove )>0)
                    
                    % point per point, remove the index to the triangle list
                    outMesh.points(pointsToRemove(1),:)=[];
                    outMesh.triangles(outMesh.triangles>pointsToRemove(1))=outMesh.triangles(outMesh.triangles>pointsToRemove(1))-1;
                    nPointsRemoved = nPointsRemoved+1;
                    pointsToRemove  = pointsToRemove-1;
                    pointsToRemove(1)=[];
                    
                end
                
                % update attributes
                for i=1:numel(obj.attributes)
                    if (numel(obj.attributes(i).attribute_array)==outMesh.npoints)
                        outMesh.attributes(i).attribute_array(pointsToRemove)=[];
                    elseif (numel(obj.attributes(i).attribute_array)==outMesh.ntriangles)
                        if max(outMesh.triangles(:))>=numel(outMesh.points)
                            %outMesh.attributes(i).attribute_array(trianglesRemoved)=[];
                            disp('MeshType: this has to be implemented')
                        end
                    end
                end
                
                outMesh.npoints = size(outMesh.points,1);
                outMesh.ntriangles = size(outMesh.triangles,1);
                
            elseif (strcmp('noholes',mode))
                
                disp(['Mode ' mode ' is not implemented yet']) ;
            else
                
                disp(['ERROR: ' mode ' is not a valid option']);
            end
            
            
            
            
        end
        
        function invertOrientation(obj)
            % This function inverts the orientation of all triangles in a
            % mesh. This might be useful for example if the mesh is negative
            % oriented and one wants it positive oriented and vv.
            if size(obj.triangles,2)~=3
                disp('MeshType::invertOrientation -- ERROR: this method only works on triangular meshes' );
                return;
            end
            obj.triangles = obj.triangles(:,[1 3 2]);
            
        end
        
        %%
        function n= GetTriangleNeighbours(obj,i)
            % n gives the neighbours triangles
            vertices = obj.triangles(i,:);
            [r1,~]= find(obj.triangles==vertices(1));
            [r2,~]= find(obj.triangles==vertices(2));
            [r3,~]= find(obj.triangles==vertices(3));
            
            isfound = zeros(obj.ntriangles,3);
            isfound(r1,1)=1;
            isfound(r2,2)=1;
            isfound(r3,3)=1;
            
            n  = find(sum(isfound,2)==2);
        end
        
        %%
        function  imposeConsistentOrientation(obj)
            % This funtion imposes the same orientation to all triangloes
            % in a mesh. This orientation is the one of the first triangle,
            % thus arbitrary, and not necessarily positive.
            
            is_oriented = false(obj.ntriangles,1);
            is_oriented(1)=true;
            %set(0,'RecursionLimit',obj.ntriangles);
            
            not_oriented = 2:obj.ntriangles;
            already_oriented_fathers = [];
            while numel(not_oriented)
                oriented = find(is_oriented)';
                oriented = setdiff(oriented,already_oriented_fathers);
                if ~numel(oriented)
                    break;
                end
                for i = oriented
                    for j=obj.triangles(i,1)
                        is_oriented = obj.orientNeighbourTriangles(j,is_oriented);
                    end
                    already_oriented_fathers = [already_oriented_fathers i];
                end
                oriented = find(is_oriented);
                not_oriented = setdiff(not_oriented,oriented);
                disp([ 'There remain ' num2str(numel(not_oriented)) ])
            end
            
            
        end
        
        function  imposeConsistentOrientationFlatMesh(obj)
            % This funtion imposes the same orientation to all triangles
            % in a mesh. This orientation is the one of the first triangle,
            % thus arbitrary, and not necessarily positive.
            % The underlying work of this method is to reorient triangles
            % with negative area.
            
            side1 = obj.points(obj.triangles(:,1),:);
            side2 = obj.points(obj.triangles(:,2),:);
            side3 = obj.points(obj.triangles(:,3),:);
            
            vector1 = side2-side1;
            vector2 = side3-side1;
            
            cp = cross(vector1,vector2);
            
            %dt = dot(vector1,vector2,2);
            %area = sum(cp.^2,2);
            %ang = acos(dt);
            
            not_oriented = find(cp(:,3)<0);
            obj.triangles(not_oriented,:) = obj.triangles(not_oriented,[1 3 2]);
        end
        
        function is_oriented = orientNeighbourTriangles(obj,i,is_oriented)
            % assuming the triangle i is oriented, set the same orientation
            % in its neighbours and propagate
            
            triangle_neighbours = obj.GetTriangleNeighbours(i); % triangles that share at least one edge
            
            neighbours_not_oriented = ~is_oriented(triangle_neighbours);
            remaining_neghbours = triangle_neighbours(neighbours_not_oriented);
            
            if ~numel(remaining_neghbours )
                return;
            end
            
            vertices = obj.triangles(i,[1 2 3 1 2]);
            for r = remaining_neghbours(:)'
                if is_oriented(r)
                    continue;
                end
                vertices_r = obj.triangles(r,[1 2 3]);
                for common_vertex=1:2
                    idx = find(vertices==vertices_r(common_vertex),1,'first');
                    if ~numel(idx)
                        continue;
                    end
                    if vertices_r(common_vertex)==vertices(idx+1)
                        %  realign it
                        obj.triangles(r,:)=obj.triangles(r,[1 3 2]);
                    end
                    break;
                end
                is_oriented(r)=true;
                % is_oriented = obj.orientNeighbourTriangles(r,is_oriented);
                
            end
            
        end
        
        %%
        function  nnonpos = imposePositiveOrientation(obj)
            % This function reorganises topology do that every facet has a
            % positive orientation. This orientation is required by many
            % algorithms and is in general desirable.
            % Functionning is only guaranteed if the gravity center of the
            % shape lays inside the mesh
            
            centre = mean(obj.points,1);
            normals = obj.GetNormalAtFaces(1:obj.ntriangles);
            inwardsVector = ones(obj.ntriangles,1)*centre   - obj.points(obj.triangles(:,1),:);
            
            isNotPositiveOriented = dot(normals,inwardsVector,2)>0;
            % take the ones not positively oriented and reorder (it suffice to swa cols 2 and 3)
            nnonpos = numel(nonzeros(isNotPositiveOriented));
            obj.triangles(isNotPositiveOriented,[2 3]) = obj.triangles(isNotPositiveOriented,[3 2]);
        end
        
        function flag = isInside(obj,points)
            % Returns a boolean indicating whether tha input point(s) is inside or outside the closed mesh
            % points should be a N x 3 matrix.
            tree = AABB_tree(obj);
            tree.BuildTree();
            
            obj.ComputeBounds();
            n_points = size(points,1);
            flag = zeros(n_points,1);
            % do a quick bounds check
            values_out = find( (points(:,1) < obj.bounds(1)) | (points(:,1) > obj.bounds(2)) |...
                (points(:,2) < obj.bounds(3)) | (points(:,2) > obj.bounds(4)) |...
                (points(:,3) < obj.bounds(5)) | (points(:,3) > obj.bounds(6)) );
            values_in = setdiff(1:n_points,values_out);
            
            % K = 2*norm(obj.bounds([2 4 6])-obj.bounds([1 3 5]));
            direction = [0 -1 0]; % any random direction would do
            
            tree.computedIntersections = zeros(numel(values_in),1);
            tree.intersect_ray(tree.rootNode,points(values_in,:), direction);
            nintersections1 = tree.computedIntersections==1; % only inside points intersect once
            
            
            flag(values_in)=nintersections1(:)';
            
            
        end
        
        
        
        function [p indicesClosestTriangle] = findClosestSurfacePoint(obj, pointIn)
            % Finds the closest point on the surface to a point (or points)
            % which in general is not a vertex.
            %1. Find closest triangle
            
            % THIS HAS TO BE REDONE
            disp('WARNING: findClosestSurfacePoint is inaccurate')
            pointsOfTriangles = cat(3,obj.points(obj.triangles(:,1),:),obj.points(obj.triangles(:,2),:),obj.points(obj.triangles(:,3),:));
            medicenters = mean(pointsOfTriangles,3);
            
            errorx = medicenters(:,1)*ones(1,size(pointIn,1)) - ones(obj.ntriangles,1) * pointIn(:,1)';
            errory = medicenters(:,2)*ones(1,size(pointIn,1)) - ones(obj.ntriangles,1) * pointIn(:,2)';
            errorz = medicenters(:,3)*ones(1,size(pointIn,1)) - ones(obj.ntriangles,1) * pointIn(:,3)';
            
            totalError = errorx.^2+errory.^2+errorz.^2;
            
            [dummy indicesClosestTriangle] = min(totalError,[],1);
            medicenters = medicenters(indicesClosestTriangle,:);
            
            
            [indicesClosestTriangle pointInTriangle]= obj.findClosestTriangle(pointIn);
            % get normals to faces
            normals = obj.GetNormalAtFaces(indicesClosestTriangle);
            p = pointIn + normals.*(dot(pointInTriangle-pointIn,normals,2) * ones(1,3));
            
        end
        
        function [indicesClosestTriangle medicenters]= findClosestTriangle(obj, pointIn)
            % Find the closest triangle to a point (or points). This is
            % usefull if for example the triangularization needs to be
            % changed
            
            % compute medicenters
            pointsOfTriangles = cat(3,obj.points(obj.triangles(:,1),:),obj.points(obj.triangles(:,2),:),obj.points(obj.triangles(:,3),:));
            medicenters = mean(pointsOfTriangles,3);
            
            errorx = medicenters(:,1)*ones(1,size(pointIn,1)) - ones(obj.ntriangles,1) * pointIn(:,1)';
            errory = medicenters(:,2)*ones(1,size(pointIn,1)) - ones(obj.ntriangles,1) * pointIn(:,2)';
            errorz = medicenters(:,3)*ones(1,size(pointIn,1)) - ones(obj.ntriangles,1) * pointIn(:,3)';
            
            totalError = errorx.^2+errory.^2+errorz.^2;
            
            [dummy indicesClosestTriangle] = min(totalError,[],1);
            medicenters = medicenters(indicesClosestTriangle,:);
            
        end
        
        
        function indicesNclosestPoints = findNClosestVertices(obj, pointIn,N)
            % Finds the N closest vertices to a point (or points)
            % PointIn has to be a N x 3 array of points
            % Note that this does not necessarily give the three vertex of
            % the closest triangle!
            
            errorx = obj.points(:,1)*ones(1,size(pointIn,1)) - ones(obj.npoints,1) * pointIn(:,1)';
            errory = obj.points(:,2)*ones(1,size(pointIn,1)) - ones(obj.npoints,1) * pointIn(:,2)';
            errorz = obj.points(:,3)*ones(1,size(pointIn,1)) - ones(obj.npoints,1) * pointIn(:,3)';
            
            totalError = errorx.^2+errory.^2+errorz.^2;
            
            [totalErrorSorted IX] = sort(totalError,1);
            
            indicesNclosestPoints = IX(1:N,:);
            
        end
        
        function [p distance] = findClosestVertex(obj, pointIn)
            % Finds the closest vertex to a point (or points)
            % PointIn has to be a N x 3 array of points
            
            errorx = obj.points(:,1)*ones(1,size(pointIn,1)) - ones(obj.npoints,1) * pointIn(:,1)';
            errory = obj.points(:,2)*ones(1,size(pointIn,1)) - ones(obj.npoints,1) * pointIn(:,2)';
            errorz = obj.points(:,3)*ones(1,size(pointIn,1)) - ones(obj.npoints,1) * pointIn(:,3)';
            
            totalError = errorx.^2+errory.^2+errorz.^2;
            
            [distance p] = min(totalError);
            
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
        
        function at = addVertexAttribute(obj, array,name)
            at = AttributeType();
            at.attribute='field';
            at.name=name;
            at.numComp=1; % between 1 and 4
            at.type='float';
            at.nelements=obj.npoints; %number of vertex
            at.attribute_array = array;
            
            % ovewrite if already exists
            ind = obj.find_attribute(at.name);
            if (ind > 0)
                obj.attributes(ind)=at;
            else
                n_attributes = numel(obj.attributes);
                obj.attributes(n_attributes+1)=at;
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
                obj.triangles(obj.triangles == points_to_be_removed(i))=replacements(i);
            end
            
            % remove points and update topology
            
            while(numel(points_to_be_removed))
                obj.triangles(obj.triangles>points_to_be_removed(1))=obj.triangles(obj.triangles>points_to_be_removed(1))-1;
                points_to_be_removed(1)=[];
                points_to_be_removed=points_to_be_removed-1;
            end
            
            obj.npoints = size(obj.points,1);
            
        end
        
        
        function propagateLabel(obj, ind, varargin)
            % This function will need some morphomat operations on meshes.
            % Then we will do a reconstruction. The seed for the
            % reconstruction can be the first point inside the region. This
            % should be easy to get since edges are oriented
            
            % the neighbouring vertices. For this, only vertices which lie
            % on the positive side of a labelled edge (edge containing two labelled
            % vertices) are labelled. Positive means on the sense of the
            % cross product of the oriented edge with the face normal.
            %
            % the label is supposed to be binary!
            
            
            
            disp('ERRROR: STILL HAS TO BE IMPLEMENTED')
            center= [];
            
            for i=1:size(varargin,2)
                if (strcmp(varargin{i},'insidePoint'))
                    center=varargin{i+1};
                end
                
            end
            
            old_labelledVertices = find(obj.attributes(ind).attribute_array);
            %  figure; viewMesh(obj,'wireframeSurface','color',[0 0 1],'labelColor',2)
            %  hold on; plot3(obj.points(old_labelledVertices,1),obj.points(old_labelledVertices,2), obj.points(old_labelledVertices,3),'.g','MarkerSize',25 ); hold off
            %  hold on; text(obj.points(:,1),obj.points(:,2), obj.points(:,3),num2str((1:obj.npoints)'),'FontSize',15 ); hold off
            for vert=old_labelledVertices'
                % get labelled neighbours
                neighbours = obj.GetForwardVertexNeighbours(vert); % this is already positively oriented
                for n = 1:size(neighbours,1)
                    if obj.attributes(ind).attribute_array(neighbours(n,1))
                        obj.attributes(ind).attribute_array(neighbours(n,2))=1;
                    end
                end
                
            end
            labelledVertices = find(obj.attributes(ind).attribute_array);
            % figure; viewMesh(obj,'wireframeSurface','color',[0 0 1],'labelColor',2)
            %hold on; plot3(obj.points(labelledVertices,1),obj.points(labelledVertices,2), obj.points(labelledVertices,3),'.g','MarkerSize',25 ); hold off
        end
        
        
        %           function propagateLabel(obj, ind, varargin)
        %             % This function propagates the label indicated by ind towards
        %             % the neighbouring vertices. For this, only vertices which
        %             lie
        %             % on the positive side of a labelled edge (edge containing two labelled
        %             % vertices) are labelled. Positive means on the sense of the
        %             % cross product of the oriented edge with the face normal.
        %             %
        %             % the label is supposed to be binary!
        %
        %
        %             center= [];
        %
        %             for i=1:size(varargin,2)
        %                 if (strcmp(varargin{i},'insidePoint'))
        %                     center=varargin{i+1};
        %                 end
        %
        %             end
        %
        %             old_labelledVertices = find(obj.attributes(ind).attribute_array);
        %             figure; viewMesh(obj,'wireframeSurface','color',[0 0 1],'labelColor',2)
        %             hold on; plot3(obj.points(labelledVertices,1),obj.points(labelledVertices,2), obj.points(labelledVertices,3),'.g','MarkerSize',25 ); hold off
        %             hold on; text(obj.points(:,1),obj.points(:,2), obj.points(:,3),num2str((1:obj.npoints)'),'FontSize',15 ); hold off
        %             for vert=labelledVertices'
        %                 % get labelled neighbours
        %                 neighbours = obj.GetForwardVertexNeighbours(vert); % this is already positively oriented
        %                 for n = 1:size(neighbours,1)
        %                     if obj.attributes(ind).attribute_array(neighbours(n,1))
        %                             obj.attributes(ind).attribute_array(neighbours(n,2))=1;
        %                     end
        %                 end
        %
        %             end
        %             labelledVertices = find(obj.attributes(ind).attribute_array);
        %             figure; viewMesh(obj,'wireframeSurface','color',[0 0 1],'labelColor',2)
        %             hold on; plot3(obj.points(labelledVertices,1),obj.points(labelledVertices,2), obj.points(labelledVertices,3),'.g','MarkerSize',25 ); hold off
        %         end
        %
        
        function [isBorder,borderPaths] = detectBorders(obj)
            % This function returns wether an vertex is a border or not
            obj.removeDuplicatedPoints();
            isBorder = zeros(obj.npoints,1);
            borderPaths=[];
            pointlist = reshape(obj.triangles,obj.ntriangles*3,1);
            pointlist = unique(pointlist);
            
            
            edges = [];
            for i=pointlist(:)'
                % get all triangles in which this point appears
                [x, ~]=find(obj.triangles == i);
                if numel(x>0)
                    neighbors = obj.triangles(x,:);
                    neighbors = neighbors(:)';
                    for j=unique(neighbors)
                        if nnz(neighbors==j)<2
                            isBorder(i,:)=1;
                            if i~=j
                                [f1,~] = find(edges==j);
                                [f2,~] = find(edges==i);
                                if ~numel(intersect(f1,f2))
                                    edges = [edges; i j];
                                end
                            end
                        end
                    end
                end
                
            end
            edges = unique(edges,'rows');
            % add attribute
            at = AttributeType();
            at.attribute='field';
            at.name='border';
            at.numComp=1; % between 1 and 4
            at.type='short';
            at.nelements=obj.npoints; %number of tuples
            for i=1:at.nelements
                at.attribute_array(i,:) = isBorder(i);
            end
            ind = obj.find_attribute(at.name);
            if (ind > 0)
                obj.attributes(ind)=at;
            else
                n_attributes = numel(obj.attributes);
                obj.attributes(n_attributes+1)=at;
            end
            
            
            
            
            % get the borderPaths
            borderPoints = find(isBorder);
            npaths = 0;
            while numel(edges );
                npaths = npaths+1;
                borderPaths{npaths} = [];
                
                currentIndex=1;
                finish=false;
                while ~finish
                    if ~numel(currentIndex)
                        finish=true;
                        continue;
                    end
                    current = edges(currentIndex,:);
                    edges(currentIndex,:)=[];
                    if numel(borderPaths{npaths}) && borderPaths{npaths}(end)==current(2)
                        current = current([2 1]);
                    end
                    extremum = current(2);
                    borderPaths{npaths} = [ borderPaths{npaths} ;current];
                    
                    % find neighbours (border edges)
                    [tri,~]=find(edges==extremum);
                    if ~numel(tri)
                        
                        % turn around starting from the other end
                        borderPaths{npaths}(:) = borderPaths{npaths}(end:-1:1);
                        current = borderPaths{npaths}(end,:);
                        [tri,~]=find(edges==current(2));
                        if ~numel(tri)
                            
                            finish=true;
                            %npaths = npaths-1;
                            %borderPaths={borderPaths{1:npaths}};
                            continue;
                        end
                    end
                    % find the edge that follows
                    e = edges(tri,:);
                    % see if we have finished
                    if numel(borderPaths{npaths})>2 && current(1)~=borderPaths{npaths}(2) && nnz(e==borderPaths{npaths}(1))
                        finish=true;
                        if e(2)==borderPaths{npaths}(end) && e(1)==borderPaths{npaths}(1) || ...
                           e(1)==borderPaths{npaths}(end) && e(2)==borderPaths{npaths}(1)
                            % this is the connecting edge - just remove it
                            edges(tri,:)=[];
                        end
                    else
                        neighbour = setdiff(e(:) ,borderPaths{npaths}(:));
                        if numel(neighbour) > 1
                            disp('MeshType::ERROR: more than one neighbour');
                            return;
                        elseif numel(neighbour)<1
                            disp('MeshType::ERROR: less than one neighbour');
                            return;
                        end
                        
                        [tri1,~]=find(edges==neighbour);
                        [tri2,~]=find(edges==current(2));
                        currentIndex = intersect(tri1,tri2);
                    end
                end
                
                
            end
            %             borderPoints = find(isBorder);
            %             npaths = 0;
            %             while numel(edges );
            %                 npaths = npaths+1;
            %                 borderPaths{npaths} = [];
            %
            %                 currentIndex=1;
            %                 finish=false;
            %                 while ~finish
            %                     if ~numel(currentIndex)
            %                         finish=true;
            %                         continue;
            %                     end
            %                     current = edges(currentIndex,:);
            %                     edges(currentIndex,:)=[];
            %                     borderPaths{npaths} = [ borderPaths{npaths} ;current];
            %
            %                     % find neighbours (border edges)
            %                     [tri,~]=find(edges==current(2));
            %                     if ~numel(tri)
            %                         finish=true;
            %                         %npaths = npaths-1;
            %                         %borderPaths={borderPaths{1:npaths}};
            %                         continue;
            %                     end
            %                     % find the edge that follows
            %                     e = edges(tri,:);
            %                     % see if we have finished
            %                     if numel(borderPaths{npaths})>2 && current(1)~=borderPaths{npaths}(2) && nnz(e==borderPaths{npaths}(1))
            %                         finish=true;
            %                     else
            %                         neighbour = setdiff(e(:) ,borderPaths{npaths}(:));
            %                         if numel(neighbour) > 1
            %                             disp('MeshType::ERROR: more than one neighbour');
            %                             return;
            %                         elseif numel(neighbour)<1
            %                             disp('MeshType::ERROR: less than one neighbour');
            %                             return;
            %                         end
            %
            %                         [tri1,~]=find(edges==neighbour);
            %                         [tri2,~]=find(edges==current(2));
            %                         currentIndex = intersect(tri1,tri2);
            %                     end
            %                 end
            %
            %
            %             end
            
        end
        
    end % methods
    
    
    
    
end

function n = cyclicNext(a,b)
n = a+1;
if (n>b)
    n = 1;
end
end

