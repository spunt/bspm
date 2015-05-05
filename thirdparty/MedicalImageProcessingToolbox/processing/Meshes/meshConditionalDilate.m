function dilatedMesh = meshConditionalDilate(meshIn, seed, varargin)
% Do conditional dilate one mesh seeded by another mesh.
% If not specified otherwise, the first point-wise label will be used in
% both cases
% The dilation is restricted to the seed


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
while (~numel(labelIn))
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
while (~numel(labelSeed))
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

dilatedMesh=MeshType(meshIn);
dilatedMesh.points = meshIn.points;
dilatedMesh.triangles = meshIn.triangles;
dilatedMesh.attributes(labelIn) = AttributeType(dilatedMesh.npoints);
dilatedMesh.attributes(labelIn).attribute_array = meshIn.attributes(labelIn).attribute_array;


 
  difference = (dilatedMesh.attributes(labelIn).attribute_array-seed.attributes(labelSeed).attribute_array)'*...
      (dilatedMesh.attributes(labelIn).attribute_array-seed.attributes(labelSeed).attribute_array);
 
while(difference)
  
     old_attribute = dilatedMesh.attributes(labelIn).attribute_array;
    
    dilatedMesh = meshDilate(dilatedMesh,'labelIndex', labelIn);
    dilatedMesh.attributes(labelIn).attribute_array = ...
       min([ dilatedMesh.attributes(labelIn).attribute_array   seed.attributes(labelSeed).attribute_array;],[],2);
    
   difference =  (old_attribute-dilatedMesh.attributes(labelIn).attribute_array)'*...
       (old_attribute-dilatedMesh.attributes(labelIn).attribute_array);
   
  
    
end
       
    



end