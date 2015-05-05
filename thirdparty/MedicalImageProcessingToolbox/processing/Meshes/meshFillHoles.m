function filledMesh=meshFillHoles(meshIn, seed, varargin)
% Do morphological hole filling  of one mesh seeded by another mesh.
% If not specified otherwise, the first point-wise label will be used in
% both cases. The seed is the index of a vertex


labelIndex= [];
i=1;
th=0;
atname='filled';
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'labelIn'))
        labelIndex = varargin{i+1};
        i = i+1;
    elseif (strcmp( varargin{i} , 'th'))
        th = varargin{i+1};
        i = i+1;
    elseif (strcmp( varargin{i} , 'name'))
        atname = varargin{i+1};
        i = i+1;
    end
    i = i+1;
end

% find the first point wise labels
if ~numel(labelIndex)
    labelIndex=-1;

    for labelIndex_=1:numel(meshIn.attributes);
        if (numel(meshIn.attributes(labelIndex_).attribute_array)==meshIn.npoints)
            labelIndex=labelIndex_;
        end
    end

    if (labelIndex<0)
        disp('ERROR: there is no attribute available')
        return;
    end
end
    

% Do the stuff. 

filledMesh=MeshType(meshIn);
filledMesh.points = meshIn.points;
filledMesh.triangles = meshIn.triangles;
filledMesh.attributes(labelIndex) = AttributeType(filledMesh.npoints);
filledMesh.attributes(labelIndex).attribute_array(seed)=1;


 
  difference = (filledMesh.attributes(labelIndex).attribute_array-meshIn.attributes(labelIndex).attribute_array)'*...
      (filledMesh.attributes(labelIndex).attribute_array-meshIn.attributes(labelIndex).attribute_array);
 
while(difference)
  
     old_attribute = filledMesh.attributes(labelIndex).attribute_array;
    
    filledMesh = meshDilate(filledMesh,'labelIndex', labelIndex);
%     filledMesh.attributes(labelIn).attribute_array = ...
%        min([ dilatedMesh.attributes(labelIn).attribute_array   meshIn.attributes(labelIndex).attribute_array;],[],2);
    filledMesh.attributes(labelIndex).attribute_array(meshIn.attributes(labelIndex).attribute_array>th)=old_attribute(meshIn.attributes(labelIndex).attribute_array>th);
       
   
   difference =  (old_attribute-filledMesh.attributes(labelIndex).attribute_array)'*...
       (old_attribute-filledMesh.attributes(labelIndex).attribute_array);
   
  
    
end

filledMesh.attributes(labelIndex).name = atname;

end