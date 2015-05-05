classdef VectorImageType < ImageType
    
    properties(GetAccess = 'public', SetAccess = 'public')
        datax=[];
        datay=[];
        dataz=[];
    end
    
    methods(Access = public)
        %constructor
        function obj = VectorImageType(size,origin,spacing,orientation)
            %obj = obj@ImageType(size,origin,spacing,orientation);
            if (nargin==1)
                % argument "size" is another images
                obj = VectorImageType(size.size,size.origin,size.spacing,size.orientation);
            elseif(nargin > 0)
                obj.size = size(:);
                obj.spacing =  spacing(:);
                obj.origin = origin(:);
                obj.orientation = orientation;
                s = obj.size;
                obj.data = zeros(s(:)');
                obj.datax = zeros(s(:)');
                obj.datay = zeros(s(:)');
                obj.dataz = zeros(s(:)');
                obj.paddingValue = s*0;
                obj.D = orientation;
                for i=1:numel(obj.size)
                    obj.D(:,i)=obj.D(:,i)*obj.spacing(i);
                end
            end
        end
        
        
        function P = GetValue(obj, pos, varargin)
            % get the pixel value (vector) at a non index position
            % can use different interpolation schemes:
            % im.GetValue(pos)  returns the value using nearest neighrbor
            % interpolation
            %   pos = world coordinates of the position where the value is
            %   desired
            % im.GetValue(pos, mode) uses the following interpolation
            %   mode = 'NN'     nearest neighbor
            %   mode = 'linear'    (tri) linear interpolation
            %   mode = 'spline'    (cubic) b-spline interpolation
            % im.GetValue(pos, mode,field) returns a scalar with the field
            % that can be:
            %   field = 'data'
            %   field = 'datax'
            %   field = 'datay'
            %   field = 'dataz'
            
            P=obj.paddingValue;
            mode = 'NN'; % NN
            field = '';
            if  (size(varargin,2)>0)
                mode = varargin{1};
            end
            if  (size(varargin,2)>1)
                field= varargin{2};
            end
            
            index = obj.GetContinuousIndex(pos);
            
            round_index = round(index);
            
            % find the indexes inside the range
            c_in = zeros(1,size(index,2));
            ndims = numel(obj.size);
            for i=1:ndims
                c_in = c_in | round_index(i,:)<1 | round_index(i,:)>obj.size(i);
            end
            
            
            if (strcmp(mode,'NN'))
                round_index(:,c_in) = round_index(:,c_in)*0+1;
                
                in_1D = sub2ind(obj.size',round_index(1,:),round_index(2,:),round_index(3,:));
                if numel(field) &&   any(strcmp(properties(obj), field))
                    P = obj.(field)(in_1D);
                else
                    P = [obj.datax(in_1D);obj.datay(in_1D);obj.dataz(in_1D)];
                end
                
            elseif (strcmp(mode,'linear'))
                
                index(:,c_in) = index(:,c_in)*0+1;
                P = obj.evaluate_at_point_linear(index,field);
            elseif (strcmp(mode,'spline'))
                index(:,c_in) = index(:,c_in)*0+1;
                P = obj.evaluate_at_point_spline(index,field);
            end
            
        end
        
        function out = magnitudeImage(obj)
            
            out = ImageType(obj);
            out.data = sqrt(obj.datax.^2+obj.datay.^2+obj.dataz);
            
        end
        
        function out = extractFrame(obj,nframe)
            if numel(obj.size)~=4
                disp('WARNING: input image is not 4D. There might be problems')
            end
            out = VectorImageType(obj.size(1:end-1),obj.origin(1:end-1),obj.spacing(1:end-1),obj.orientation(1:end-1,1:end-1));
            out.data = obj.data(:,:,:,nframe);
            out.datax = obj.datax(:,:,:,nframe);
            out.datay = obj.datay(:,:,:,nframe);
            out.dataz = obj.dataz(:,:,:,nframe);
            
        end
        
        function streamlines(obj,startpoints,varargin)
            pos = obj.GetPosition(1:prod(obj.size))';
            inplane=false;
            i=1;
            vecoptions={};
            while (i <= size(varargin,2))
                if  (strcmp( varargin{i} , 'inplane'))
                    inplane= true;
                elseif  (strcmp( varargin{i} , 'options'))
                    options= varargin{i+1};
                end
                i = i+1;
            end
            
            if numel(obj.size)==2
                
                x = reshape(pos(:,1),obj.size');
                y = reshape(pos(:,2), obj.size');
                u = reshape(obj.datax(:),obj.size');
                v = reshape(obj.datay(:),obj.size');
                
                XY_ = stream2(y,x,v,u,startpoints(:,2),startpoints(:,1),[0.2 10000]);
                id=1;
                for jj = 1:numel(XY_)
                    if numel(XY_{jj})
                        XY{id}=XY_{jj}(:,[2 1]);
                        id = id+1;
                    end
                end
                
                hlines = streamline(XY);
                set(hlines,options{1:end});
                
            elseif numel(obj.size)==3
                if inplane
                    idx =find(obj.size==1);
                    through_planeDirection = obj.orientation(:,idx);
                    tp_component = dot([obj.datax(:) obj.datay(:) obj.dataz(:)],ones(numel(obj.datax),1)*through_planeDirection',2);
                    tp_component(tp_component~=tp_component)=0;
                    dx = obj.datax(:)-tp_component*through_planeDirection(1);
                    dy = obj.datay(:)-tp_component*through_planeDirection(2);
                    dz = obj.dataz(:)-tp_component*through_planeDirection(3);
                    
                    X = reshape(pos(:,1),obj.size');
                    Y = reshape(pos(:,2),obj.size');
                    Z = reshape(pos(:,3),obj.size');
                    U = reshape(dx,obj.size');
                    V = reshape(dy,obj.size');
                    W = reshape(dz,obj.size');
                    
                    streamline(X,Y,Z,U,V,W,startpoints(:,1),startpoints(:,2),startpoints(:,3),vecoptions{1:end});
                else
                    streamline(pos(:,1),pos(:,2),pos(:,3),obj.datax(:),obj.datay(:),obj.dataz(:),startpoints(:,1),startpoints(:,2),startpoints(:,3),vecoptions{1:end});
                end
            else
                disp('streamlines method is only available for 2D or 3D images')
            end
        end
        
        function quiver(obj,varargin)
            pos = obj.GetPosition(1:prod(obj.size))';
            inplane=false;
            i=1;
            vecoptions={};
            use_jet=false;
            while (i <= size(varargin,2))
                if  (strcmp( varargin{i} , 'inplane'))
                    inplane= true;
                elseif  (strcmp( varargin{i} , 'quiver_options'))
                    vecoptions= varargin{i+1};
                elseif  (strcmp( varargin{i} , 'use_jet'))
                    use_jet= true;
                end
                i = i+1;
            end
            
            if numel(obj.size)==2
                quiver(pos(:,1),pos(:,2),obj.datax(:),obj.datay(:),vecoptions{1:end});
            elseif numel(obj.size)==3
                if inplane
                    idx =find(obj.size==1);
                    through_planeDirection = obj.orientation(:,idx);
                    tp_component = dot([obj.datax(:) obj.datay(:) obj.dataz(:)],ones(numel(obj.datax),1)*through_planeDirection',2);
                    tp_component(tp_component~=tp_component)=0;
                    dx = obj.datax(:)-tp_component*through_planeDirection(1);
                    dy = obj.datay(:)-tp_component*through_planeDirection(2);
                    dz = obj.dataz(:)-tp_component*through_planeDirection(3);
                    quiver3(pos(:,1),pos(:,2),pos(:,3),dx,dy,dz,vecoptions{1:end});
                else
                    
                    if use_jet
                        arrow_colors =jet(32);
                        obj.data = sqrt(obj.datax.^2+obj.datay.^2+obj.dataz.^2);
                        M = max(obj.data(:));
                        m = min(obj.data(:));
                        inter = (M-m)/size(arrow_colors,1);
                        hold on;
                        for i=1:size(arrow_colors,1)
                            idx = find((obj.data(:)< m+inter*i & obj.data(:)>= m+inter*(i-1)) & abs(obj.data(:))>0);
                            reoriented_vectors = (obj.orientation * [obj.datax(idx) obj.datay(idx) obj.dataz(idx)]')';
                            quiver3(pos(idx,1),pos(idx,2),pos(idx,3),reoriented_vectors(:,1),reoriented_vectors(:,2),reoriented_vectors(:,3),vecoptions{1:end},'Color',arrow_colors(i,:));
                        end
                    else
                        quiver3(pos(:,1),pos(:,2),pos(:,3),obj.datax(:),obj.datay(:),obj.dataz(:),vecoptions{1:end});
                    end
                end
            else
                disp('Quiver method is only available for 2D or 3D images')
            end
        end
        
        function show(obj,varargin)
            pos = obj.GetPosition(1:prod(obj.size))';
            i=1;
            component='data';
            while (i <= size(varargin,2))
                if  (strcmp( varargin{i} , 'component'))
                    component = varargin{i+1};
                    i = i+1;
                end
                i = i+1;
            end
            
            
            if numel(obj.size)==2
                
            elseif numel(obj.size)==3
                
                if ndims(obj.datax)==2 % this is a slice of a 3D image
                    %bounds = [min(pos) max(pos)];
                    X = reshape(pos(:,1),size(obj.datax));
                    Y = reshape(pos(:,2),size(obj.datax));
                    Z = reshape(pos(:,3),size(obj.datax));
                    if strcmp(component,'data')
                        surf(X, Y, Z, sqrt(obj.datax.^2+obj.datay.^2+obj.dataz.^2), 'FaceColor', 'texturemap','EdgeColor','none');
                    else
                        surf(X, Y, Z, obj.(component), 'FaceColor', 'texturemap','EdgeColor','none');
                    end
                else
                    disp('show method not implemented for 3D images yet')
                end
            else
                disp('show method not implemented for 3D images yet')
            end
        end
    end
    
    
    methods(Access = private)
        
        function  value = evaluate_at_point_linear(obj,continuous_index,varargin)
            % continuous_index is K x N
            field = '';
            if  (size(varargin,2)>0)
                field = varargin{1};
            end
            
            
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index(' num2str(i) ',:)'];
            end
            
            if numel(field) &&  any(strcmp(properties(obj), field))
                value = eval(['interpn(obj.data' str ',''linear'')']);
            else
                valuex = eval(['interpn(obj.datax' str ',''linear'')']);
                valuey = eval(['interpn(obj.datay' str ',''linear'')']);
                valuez = eval(['interpn(obj.dataz' str ',''linear'')']);
                
                
                value = [valuex; valuey; valuez];
            end
            
            
        end
        
        
        
        function  value = evaluate_at_point_spline(obj,continuous_index,varargin)
            % continuous_index is 3 x N
            
            field = '';
            if  (size(varargin,2)>0)
                field = varargin{1};
            end
            
            
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index(' num2str(i) ',:)'];
            end
            
            if numel(field) &&  any(strcmp(properties(obj), field))
                value = eval(['interpn(obj.data' str ',''cubic'')']);
            else
                valuex = eval(['interpn(obj.datax' str ',''cubic'')']);
                valuey = eval(['interpn(obj.datay' str ',''cubic'')']);
                valuez = eval(['interpn(obj.dataz' str ',''cubic'')']);
                
                
                value = [valuex; valuey; valuez];
            end
            
            
        end
        
    end
end