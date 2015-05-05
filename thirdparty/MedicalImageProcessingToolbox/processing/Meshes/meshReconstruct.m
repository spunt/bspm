function meshReconstruct(meshIn, seed, varargin)
% Do morphological reconstruction of one mesh seeded by another mesh.
% If not specified otherwise, the first point-wise label will be used in
% both cases
% TODO!


labelSeed = [];
labelIn= [];
i=1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'labelIn'))
        labelIn = varargin{i+2};
        i = i+1;
    elseif(strcmp( varargin{i} , 'labelSeed'))
        labelSeed= varargin{i+1};
        i = i+1;
    end
    i = i+1;
end

% find the first point wise labels
labelIndex_=1;
while (~numel(labelSeed))
    if (numel(meshIn.attributes(labelIndex_).attribute_array)==meshIn.npoints)
        labelIn=labelIndex_;
    else
        labelIndex_=labelIndex_+1;
    end
    if (labelIndex_>numel(meshIn.attributes))
        disp('ERROR: there is no attribute available')
        return
    end
end

labelIndex_=1;
while (~numel(labelIn))
    if (numel(seed.attributes(labelIndex_).attribute_array)==seed.npoints)
        labelSeed=labelIndex_;
    else
        labelIndex_=labelIndex_+1;
    end
    if (labelIndex_>numel(seed.attributes))
        disp('ERROR: there is no attribute available')
        return
    end
end

% Do the stuff. 

meshOut=MeshType(meshIn);
meshOut.points = meshIn.points;
meshOut.triangles = meshIn.triangles;

 at = AttributeType();
 at.attribute='field';
 at.name='eroded';
 at.numComp=1; % between 1 and 4
 at.type='float';
 at.nelements=meshIn.npoints; %number of vertex
 at.attribute_array = zeros(meshIn.npoints,1);

for i=1:meshIn.npoints
       % get neighbourgs
       nei = meshOut.GetVertexNeighbours(i);
       
       newVal = min([ meshIn.attributes(labelIndex).attribute_array(i) ;meshIn.attributes(labelIndex).attribute_array(nei)]);
       at.attribute_array(nei)=newVal;
       at.attribute_array(i)=newVal;

end

meshOut.attributes(numel(meshOut.attributes)+1)=at;

end