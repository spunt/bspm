function out = resample_lines_mesh(mesh,beamsource,factor, varargin)
% This function assumes that the mesh is provided as scanlines. resamples
% each scanline

interp='NN';
doblurring=false;
MAX_CHUNK_SIZE = 50;
s=[-1 -1 -1]';
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'interpolation'))
        interp=varargin{i+1};
    elseif (strcmp(varargin{i},'blur'))
        doblurring = true;
        blur_neighourhood=varargin{i+1};
        blur_sigma=varargin{i+2};
    elseif (strcmp(varargin{i},'size'))
        s=varargin{i+1};
    elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end
    
end
%----------------------------



m = MeshType(mesh);

% create point matrix

if s(1)==-1
    disp('resample_lines_mesh::Error, the size has to be specified')
    out=[];
    return;
end




if doblurring
            
            blur_neighourhood= floor(blur_neighourhood/2)*2+1;
            
            kernel = ones(1,blur_neighourhood);
            kernel(1,:)=gausswin(blur_neighourhood,2/blur_sigma);
            kernel = kernel/sum(kernel(:)); 
            
            
            
            im_.data =convn(ref.data,kernel,'same');
end


% -----------------------------Interpolate ----------------
% create the grid of the reference image
ndims = numel(ref.size);

NCHUNKS = ceil(ref.size/MAX_CHUNK_SIZE);

if (ndims==4)
    
    if (isa(im_,'VectorImageType') || isfield(im_,'datax'))
        out = VectorImageType(ref);
    else
        out = ImageType(ref);
    end
    
    chunked_size = ceil(out.size./NCHUNKS)';
    
    [ix, iy, iz, it]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1,0:NCHUNKS(4)-1);
    intervals = [ix(:) iy(:) iz(:) it(:)];
    clear ix iy iz it;
    for i=1:size(intervals,1)
        ranges([1 3 5 7]) = intervals(i,:).*chunked_size+1;
        ranges([2 4 6 8]) = min([(intervals(i,:)+[1 1 1 1]).*chunked_size ; out.size']);
        
        ranges_size = ranges([2 4 6 8])-ranges([1 3 5 7])+[1 1 1 1];
        % generate all the indexes of the target image
        [x, y, z, t] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8));
        
        positions = ref.GetPosition([x(:) y(:) z(:) t(:)]');
        clear x y z t;
        datas = im_.GetValue(positions,interp);
        
        if isa(im_,'VectorImageType') || any(strcmp(properties(im_), 'datax'))
            out.datax(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = reshape(datas(1,:),ranges_size);
            out.datay(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = reshape(datas(2,:),ranges_size);
            out.dataz(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = reshape(datas(3,:),ranges_size);
            
            data_ = im_.GetValue(positions,interp,'data');
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = ...
                reshape(data_,ranges_size);
            clear data positions;
        else
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = ...
                reshape(datas,ranges_size);
        end
        clear datas;
        
        
    end
    
elseif (ndims==3)
    %[X Y Z]=ndgrid(1:ref.size(1),1:ref.size(2),1:ref.size(3));
    % Retrieve the world coordinates on the reference image
    
    if (isa(im_,'VectorImageType') || isfield(im_,'datax'))
        out = VectorImageType(ref);
    else
        out = ImageType(ref);
    end
    
    chunked_size = ceil(out.size./NCHUNKS)';
    
    [ix, iy, iz]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1);
    intervals = [ix(:) iy(:) iz(:)];
    clear ix iy iz;
    for i=1:size(intervals,1)
        ranges([1 3 5]) = intervals(i,:).*chunked_size+1;
        ranges([2 4 6]) = min([(intervals(i,:)+[1 1 1]).*chunked_size ; out.size']);
        
        ranges_size = ranges([2 4 6])-ranges([1 3 5])+[1 1 1];
        % generate all the indexes of the target image
        [x, y, z] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6));
        
        positions = ref.GetPosition([x(:) y(:) z(:)]');
        clear x y z;
        datas = im_.GetValue(positions,interp);
        
        if isa(im_,'VectorImageType') || any(strcmp(properties(im_), 'datax'))
            out.datax(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(datas(1,:),ranges_size);
            out.datay(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(datas(2,:),ranges_size);
            out.dataz(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(datas(3,:),ranges_size);
            
            data_ = im_.GetValue(positions,interp,'data');
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = ...
                reshape(data_,ranges_size);
            clear data positions;
        else
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = ...
                reshape(datas,ranges_size);
        end
        clear datas;
        
        
    end
elseif(ndims==2)
    [X, Y]=ndgrid(1:ref.size(1),1:ref.size(2));
    % Retrieve the world coordinates on the reference image
    positions = ref.GetPosition([X(:) Y(:)]');
    
    if isa(im_,'VectorImageType') ||  any(strcmp(properties(obj), 'datax'))
        out = VectorImageType(ref);
        datas = im_.GetValue(positions,interp);
        out.datax = reshape(datas(1,:),ref.size(1),ref.size(2));
        out.datay = reshape(datas(2,:),ref.size(1),ref.size(2));
        % out.dataz = reshape(datas(3,:),ref.size(1),ref.size(2));
    else
        out = ImageType(ref);
        out.data = reshape(im_.GetValue(positions,interp),ref.size(1),ref.size(2));
    end
    
end
out.data(out.data~=out.data)=0;

end