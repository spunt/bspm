function [roinormal, roiorigin, points2D, points3D] =  read_roi(roi_file)
% [roinormal, roiorigin, points2D, points3D] =  read_roi(roi_file)
%roi_file must be an xml file
% points are wc in 2D of the roi vertices. This space can be achieved by
% the matrix defined by origin and normal

roi_file = char(roi_file);

idx = strfind(roi_file,'xml');
roiFileTxT = [ roi_file(1:idx-1) 'txt'];
planeData =dlmread(roiFileTxT);

roi.normal = planeData(1,:)';
roi.origin = planeData(2,:)';

roinormal = roi.normal;
roiorigin = roi.origin;

[x y]=vtkMathPerpendiculars(roi.normal,pi/2);

M=eye(4);
M(1:3,1:3) = [x(:) y(:) roi.normal(:)];
M(1:3,4) = roi.origin;

% --------------------------------

xDoc= xmlread(roi_file);

children = parseChildNodes(xDoc);
npoints=str2num(children.Children(2).Children(2).Attributes(2).Value);
points_string = children.Children(2).Children(2).Children(6).Children(2).Children.Data;
[token remain]=strtok(points_string,' ');
for i=1:npoints
    for j=1:3
    [token remain]=strtok(remain,' ');
    points2D{i,j}=str2num(token);
    end
end
points2D = cell2mat(points2D)';
points2D(3,:)=0.0;

points3D = M * [points2D ; ones(1,size(points2D,2))];

end


function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end
end

% ----- Subfunction MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end
end

% ----- Subfunction PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end
end
