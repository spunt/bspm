function imageOut = meshToImage(meshIn, field,  varargin)
% converts a mesh into an image (rasterization), on the grid defined by
% the spacing
%
%   Options:
%   'thickness' , n     ->  n is the number of voxel sizes at each size of
%   the mesh that will be painted (by default, n=1/2)

dbg=false;
i=1;
spacing = [1 1 1]';
MAX_CHUNK_SIZE = 4000;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'spacing'))
        spacing = varargin{i+1};
        i = i+1;
    elseif(strcmp( varargin{i} , 'debug'))
        dbg= true;
    elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end
    i = i+1;
end


%find field
labelIndex=[];
labelIndex_=1;
while (~numel(labelIndex))
    if strcmp(meshIn.attributes(labelIndex_).name,field)
        labelIndex=labelIndex_;
    else
        labelIndex_=labelIndex_+1;
    end
    if (labelIndex_>numel(meshIn.attributes))
        disp('ERROR: there is no attribute available')
        return
    end
end
% labelIndex(1) contains the first occurence of the desired field.

% automatically find image parameters
bounds([1 3 5])=min(meshIn.points);
bounds([2 4 6])=max(meshIn.points);
bounds = bounds(:);
s = ceil((bounds([2 4 6])-bounds([1 3 5])+0.5)./spacing);
s(3)=max(s(3),1);
isvector = false;

imageOut = ImageType(s,bounds([1 3 5]),spacing,eye(3));
if size(meshIn.attributes(labelIndex(1)).attribute_array,1)>1 && size(meshIn.attributes(labelIndex(1)).attribute_array,2)>1
    % this is a vector image
    imageOut = VectorImageType(s,bounds([1 3 5]),spacing,eye(3));
    isvector = true;
end



npoints = prod(s);

nintervals = ceil(npoints/MAX_CHUNK_SIZE);
all_distances = [];
for i=1:nintervals
    
    i1 = (i-1)*MAX_CHUNK_SIZE+1;
    i2 = min([i*MAX_CHUNK_SIZE npoints]);
    
    [x, y, z]=ind2sub(s',i1:i2);
    positions = imageOut.GetPosition([x(:) y(:) z(:)]');
    [indices, distances]= meshIn.findClosestVertex(positions');
    all_distances = [all_distances distances];
    if isvector
        imageOut.datax(i1:i2) = meshIn.attributes(labelIndex(1)).attribute_array(indices,1);
        imageOut.datay(i1:i2) = meshIn.attributes(labelIndex(1)).attribute_array(indices,2);
        imageOut.dataz(i1:i2) = meshIn.attributes(labelIndex(1)).attribute_array(indices,3);
    else
        imageOut.data(i1:i2) = meshIn.attributes(labelIndex(1)).attribute_array(indices);
    end
end




dist_th = 0.5*sqrt(sum(spacing.^2));

if isvector
    imageOut.datax(all_distances>=dist_th)=0;
    imageOut.datay(all_distances>=dist_th)=0;
    imageOut.dataz(all_distances>=dist_th)=0;
else
    imageOut.data(all_distances>=dist_th)=0;
end




end