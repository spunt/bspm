classdef ImageType < handle
    % This class is compliant with itk::Image class
    % by Alberto Gomez, 2011
    %
    % Update: add Support for passing several points, for speed (start 18/oct/2011)
    
    properties(GetAccess = 'public', SetAccess = 'public')
        data=[];
        size=[];
        origin=[]; % first voxel is indexed by [1 1 1] !!!
        spacing= [];
        orientation = []; % orientation matrix (3 x 3)
        D = [];
        paddingValue = 0;
        MAX_CHUNK_SIZE = 150; % for internal operations, to preserve memory, maximum operation block will be of 100 x 100 x 100
    end
        
    
    methods(Access = public)
        %constructor
        function obj = ImageType(size,origin,spacing,orientation)
            if (nargin==1)
                % argument "size" is another images
                obj = ImageType(size.size,size.origin,size.spacing,size.orientation);
                
            elseif (nargin==2) &&  strcmp(origin,'copy')
                obj = ImageType(size.size,size.origin,size.spacing,size.orientation);
                obj.data = size.data;
                
            elseif(nargin > 0)
                obj.size = size(:);
                obj.spacing =  spacing(:);
                obj.origin = origin(:);
                obj.orientation = orientation;
                s = obj.size;
                if numel(size)==1
                    obj.data = zeros(s(:)',1);
                    obj.size = [size(:) ; 1];
                else
                    obj.data = zeros(s(:)');
                end
                
                obj.D = orientation;
                
                for i=1:numel(obj.spacing)
                    obj.D(:,i)=obj.D(:,i)*obj.spacing(i);
                end
                
            end
            
            
        end
        
        function [P I] = GetBounds(obj,th,addBorder)
            % get the limits in x, y and z, in wc and in index
            %addborder is a boolean
            if (nargin==1)
                P = zeros(2*ndims(obj.data),1);
                % find the cornes
                nn = numel(obj.size);
                
                
                str1 = '';
                str2 = '';
                str3 = '';
                intervals=[];
                for i = 1:nn
                    str1 = [str1 '[1 obj.size('  num2str(i) ')],'];
                    str2 = [str2 'i' num2str(i) ' '];
                    str3 = [str3 'i' num2str(i) '(:) '];
                end
                str1=str1(1:end-1); % remove last comma
                eval([ '[' str2 ']= ndgrid(' str1 '); intervals = [' str3 ']; clear ' str2 ';']);
                intervals = unique(intervals,'rows');
                c = zeros(nn,size(intervals,1));
                for i = 1:size(intervals,1)
                    c(:,i) = obj.GetPosition(intervals(i,:)');
                end
                
                P= reshape([min(c,[],2) max(c,[],2)]',1,[]);
                
                
                I = [-1 -1 -1 -1 -1 -1];
            elseif (nargin>1)
                % get the limits in x, y and z
                % get limit in x, y, z where abs value of intensity is over
                % the threshold th
                % project into xy;
                for i=1:numel(obj.size)
                    otherdims = setdiff(1:ndims(obj.data),i);
                    data = abs(obj.data);
                    for j=otherdims
                        data = max(data,[],j);
                    end
                    bounds(i,1)=min(find(data));
                    bounds(i,2)=max(find(data));
                end
                
                
                P = obj.GetPosition(bounds);
                order = reshape([1:numel(obj.size); (1:numel(obj.size))+numel(obj.size)],2*numel(obj.size),[]);
                P=P(order);
                I = bounds(order);
                if (nargin==3)
                    if (addBorder)
                        I = I + reshape([-ones(1,ndims(obj.data)) ; ones(1,ndims(obj.data))],2*ndims(obj.data),[]); % leave a border
                        P = P + reshape([-obj.spacing' ; obj.spacing'],2*ndims(obj.data),[]); % leave a border
                    end
                end
            end
            P = P(:);
            I=I(:);
        end
        
        
        function pos = GetPosition(obj, index)
            % If obj is a K x n image, where K is the dimension (eg 2, 3 4, ...)
            % get the position of a voxel at index "index", a K x n matrix
            % also returns a K x n matrix
            %
            
            
            nindices = size(index,2);
            if size(index,1)==1 % If index is 1D, is assumes conversion from ind to sub
                str = '';
                str2 = '';
                for i = 1:numel(obj.size)
                    str = [str 'i'  num2str(i) ','];
                    str2 = [str2 'i'  num2str(i) '(:) '];
                end
                
                eval([ '[' str(1:end-1) ']=ind2sub(obj.size'',index); index =[' str2(1:end-1) ']''; ']);
            end
            
            pos = zeros(size(index));
            
            if nindices > obj.MAX_CHUNK_SIZE^size(index,1)
                % Divide into chunks
                NCHUNKS = ceil(nindices/( obj.MAX_CHUNK_SIZE^size(index,1)));
                chunked_size = ceil(nindices/NCHUNKS);
                intervals= 0:NCHUNKS-1;
                
                for i=1:numel(intervals)
                    ranges(1) = intervals(i)*chunked_size+1;
                    ranges(2) = min([(intervals(i)+1)*chunked_size ; nindices]);
                    nranges = ranges(2)-ranges(1)+1;
                    pos(:,ranges(1):ranges(2)) = obj.origin(:) * ones(1,nranges) + obj.D *(index(:,ranges(1):ranges(2))-1);
                end
                
            else
                pos = obj.origin(:) * ones(1,nindices) + obj.D *(index-1);
            end
        end
        
        function index = GetContinuousIndex(obj, pos)
            % get the continuous index position
            % pos has to be a K x n matrix where K is the dimensionality of
            % the image
            
            index =  obj.D\(pos - obj.origin(:)*ones(1,size(pos,2)) ) + 1;
        end
        
        function P = GetPixel(obj, index)
            % get the pixel value at an index position. "index" is a 3 x 1 matrix
            P = obj.data(index(1),index(2),index(3));
        end
        
        function c_in = isOutOfRange(obj,index)
            %  P = isInRange(index)
            % index must be a 3XN array
            % returns 0 if the index is inside the image, 0 otherwise
            
            
            c_in = zeros(1,size(index,2));
            ndims = numel(obj.size);
            for i=1:ndims
                c_in = c_in | index(i,:)<1 | index(i,:)>obj.size(i);
            end
            
        end
        
        function P = GetValue(obj, pos, varargin)
            % get the pixel value at a non index position
            % can use different interpolation schemes:
            % im.GetValue(pos)  returns the value using nearest neighrbor
            % interpolation, with pos given in world coordinates
            % im.GetValue(pos, mode) uses the following interpolation
            %   mode = 'NN'     nearest neighbor
            %   mode = 'linear'    (tri) linear interpolation
            %   mode = 'spline'    (cubic) b-spline interpolation
            
            P=obj.paddingValue;
            mode = 'NN'; % NN
            ndims = numel(obj.size);
            if  (size(varargin,2)>0)
                mode = varargin{1};
            end
            
            
            
            index = obj.GetContinuousIndex(pos);
            
            
            round_index = round(index);
            
            % find the indexes inside the range
            
            c_in = obj.isOutOfRange(round_index);
            
            if numel(c_in~=0)
                
                
                if (strcmp(mode,'NN'))
                    round_index(:,c_in) = round_index(:,c_in)*0+1;
                    
                    str = '';
                    
                    for i = 1:ndims
                        str = [str ', round_index('  num2str(i) ',:)'];
                    end
                    
                    in_1D = eval(['sub2ind(obj.size''' str ')']);
                    P = obj.data(in_1D);
                    
                elseif (strcmp(mode,'linear'))
                    
                    index(:,c_in) = index(:,c_in)*0+1;
                    P = obj.evaluate_at_point_linear(index);
                elseif (strcmp(mode,'spline'))
                    index(:,c_in) = index(:,c_in)*0+1;
                    P = obj.evaluate_at_point_spline(index);
                end
                P(c_in)=0;
                
            end
            
            
            
            
        end
        
        function SetOriginToCenter(obj)
            % TODO apply orientation here!
            obj.origin = obj.orientation * ((obj.size -1).*obj.spacing/2.0 );
        end
        
        function out = extractFrame(obj,nframe)
            if numel(obj.size)~=4
                disp('WARNING: input image is not 4D. There might be problems')
            end
            out = ImageType(obj.size(1:end-1),obj.origin(1:end-1),obj.spacing(1:end-1),obj.orientation(1:end-1,1:end-1));
            out.data = obj.data(:,:,:,nframe);
            
        end
        
        
        
        function h = show(obj,varargin)
            pos = obj.GetPosition(1:prod(obj.size))';
            opacity=-1;
            i=1;
            crange = [min(obj.data(:)) max(obj.data(:))];
            while (i <= size(varargin,2))
                if  (strcmp( varargin{i} , 'opacity'))
                    opacity = varargin{i+1};
                    i = i+1;
                elseif  (strcmp( varargin{i} , 'colorrange'))
                    crange  = varargin{i+1};
                    i = i+1;
                end
                i = i+1;
            end
            
            
            if numel(obj.size)==2
                X = reshape(pos(:,1),size(obj.data));
                Y = reshape(pos(:,2),size(obj.data));
                Z = zeros(size(obj.data));
                
                data_to_plot=(double(obj.data)-min(obj.data(:)))/(max(obj.data(:))-min(obj.data(:))); % normalised [0,1]
                data_to_plot = data_to_plot*(crange(2)-crange(1))+crange(1);
                
                if opacity>=0
                    h= surf(X, Y, Z, data_to_plot, 'FaceColor', 'texturemap','EdgeColor','none','FaceAlpha',opacity,'CDataMapping','direct');
                else    
                    h = surf(X, Y, Z, data_to_plot, 'FaceColor', 'texturemap','EdgeColor','none','CDataMapping','direct');
                end
            elseif numel(obj.size)==3
                
                if ndims(obj.data)==2 % this is a slice of a 3D image
                    %bounds = [min(pos) max(pos)];
                    X = reshape(pos(:,1),size(obj.data));
                    Y = reshape(pos(:,2),size(obj.data));
                    Z = reshape(pos(:,3),size(obj.data));
                    if opacity>=0
                       h = surf(X, Y, Z, obj.data, 'FaceColor', 'texturemap','EdgeColor','none','FaceAlpha',opacity);
                        
                    else
                    h= surf(X, Y, Z, obj.data, 'FaceColor', 'texturemap','EdgeColor','none');
                    end
                else
                    
                    bds = obj.GetBounds();
                    centroid = (bds([1 3 5])+bds([2 4 6]))/2;
                    n1 = resliceImage(obj,'plane',[1 0 0]',centroid);
                    n2 = resliceImage(obj,'plane',[0 1 0]',centroid);
                    n3 = resliceImage(obj,'plane',[0 0 1]',centroid);
                    
                    hold on; 
                    h(1)=n1.show();
                    h(2)=n2.show();
                    h(3)=n3.show();
                    hold off;
                    
                    
                end
            else
                disp('show method not implemented for 4D images yet')
            end
        end
        
    end
    
    methods(Access = private)
        function  value = evaluate_at_point_linear(obj,continuous_index)
            % continuous_index is 3 x N
            
            ind = find(obj.size<2);
            
            
            str = '';
            for i = 1:size(continuous_index,1)
                str = [str ', continuous_index('  num2str(i) ',:)'];
            end
            
            
            
            if ~numel(ind)
                % there are no singleton dimensions
                value = eval(['interpn(double(obj.data)' str ',''linear'')']);
            else
                if numel(ind)==1
                    data2 = cat(ind,obj.data,obj.data,obj.data);
                    continuous_index(ind,:)=continuous_index(ind,:)+1;
                    value = eval(['interpn(double(data2)' str ',''linear'')']);
                end
            end
           
            
        end
        
        function  value = evaluate_at_point_spline(obj,continuous_index)
            % continuous_index is K x N
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index('  num2str(i) ',:)'];
            end
            
            value = eval(['interpn(double(obj.data)' str ',''cubic'')']);
        end
        
    end % private methods
end

