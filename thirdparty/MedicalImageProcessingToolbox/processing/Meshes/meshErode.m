function meshOut = meshErode(meshIn, varargin)
% Permorms mathematicall morphology on the surface.
% Connectivity is given by topology.
% If not specified otherwise, first point-wise label is used as binary
% image
%
% Options:  'label', labelName
% Options:  'labelIndex', labelIndex

labelIndex = [];
i=1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'label'))
        labelIndex = meshIn.findAttribute(varargin{i+2})
        i = i+1;
    elseif(strcmp( varargin{i} , 'labelIndex'))
        labelIndex= varargin{i+1};
        i = i+1;
    end
    i = i+1;
end
labelIndex_=1;
while (~numel(labelIndex))
    if (numel(meshIn.attributes(labelIndex_).attribute_array)==meshIn.npoints)
        labelIndex=labelIndex_;
    else
        labelIndex_=labelIndex_+1;
    end
    if (labelIndex_>numel(meshIn.attributes))
        disp('ERROR: there is no attribute available')
        return
    end
end

meshOut=MeshType(meshIn);
meshOut.points = meshIn.points;
meshOut.triangles = meshIn.triangles;

 at = AttributeType();
 at.attribute='field';
 at.name='eroded';
 at.numComp=1; % between 1 and 4
 at.type='float';
 at.nelements=meshIn.npoints; %number of vertex
 at.attribute_array = ones(meshIn.npoints,1);

for i=1:meshIn.npoints
       % get neighbourgs
       nei = meshOut.GetVertexNeighbours(i);
       newVal = min( meshIn.attributes(labelIndex).attribute_array([ i; nei]));
       at.attribute_array( i)=newVal;
end

meshOut.attributes(labelIndex)=at;


end