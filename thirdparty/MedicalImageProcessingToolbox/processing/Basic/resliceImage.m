function [slice,slice2D] = resliceImage(im,varargin)
% out = resliceImage(im,mat)
% out = resampleImage(im,mat, options)
%
%   Options:
%       'interpolation',interp
%       'spacing',spacing (of the resulting 2D slice)
%       'mat',mat uses a matrix for slicing (4x4)
%       'plane',planenormal,planepoint
%       'roi',roifile.xml
%       'maskwithroi'       if selected, will remove everithing outside the
%       roi
%
%   NOTE: select only one of 'mat', 'plane' or 'roi'!
%
%   slice will be the slice in 3D, slice2D is a 2D image. If it is a vector
%   image, slice2D contains only the inplane velocities
%   Works for ImageType and VectorImageType


thicknes = 1;
interpolation='NN';
spacing = [-1 -1 -1];
usematrix = false;
useplane=false;
useroi=false;
doPlot=false;
maskWithRoi=false;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'interpolation'))
        interpolation=varargin{i+1};
    elseif (strcmp(varargin{i},'spacing'))
        spacing=varargin{i+1};
    elseif (strcmp(varargin{i},'mat'))
        M=varargin{i+1};
        usematrix=true;
    elseif (strcmp(varargin{i},'thicknes'))
        thicknes=varargin{i+1};
    elseif (strcmp(varargin{i},'plane'))
        slicing_plane.n=varargin{i+1};
        slicing_plane.p=varargin{i+2};
        useplane=true;
    elseif (strcmp(varargin{i},'roi'))
        roi_file=varargin{i+1};
        useroi=true;
    elseif (strcmp(varargin{i},'maskWithRoi'))
        maskWithRoi=true;
    end
    
end
%----------------------------

if (usematrix)
    slicing_plane.n = M(1:3,3);
    slicing_plane.p = M(1:3,4);
elseif (useplane)
    [xN, yN ]=vtkMathPerpendiculars(slicing_plane.n,pi/2);
    M = [xN yN slicing_plane.n slicing_plane.p ; [0 0 0 1] ];
elseif(useroi)
    [roinormal, roiorigin, roipoints] =  read_roi(roi_file);
    
    slicing_plane.n = roinormal(:);
    slicing_plane.p = roiorigin(:);
    
    [x, y]=vtkMathPerpendiculars( slicing_plane.n,pi/2);
    
    M=eye(4);
    M(1:3,1:3) = [x(:) y(:)  slicing_plane.n];
    M(1:3,4) =slicing_plane.p;
    
end

%% find the bounding box of the output slice, in 2D :
% get bound planes of image
bounds = im.GetBounds();

[x, y, z]=ndgrid(0:1,0:1,0:1);

boundsIndices = [x(:) y(:) z(:)]'+ [1 3 5]'*ones(1,8);
pointsBounds = zeros(3,8);
for i=1:8
    pointsBounds(:,i)    = bounds(boundsIndices(:,i) );
end

pointsFaces(:,1:3)=pointsBounds(:,1)*ones(1,3);
pointsFaces(:,4:6)=pointsBounds(:,end)*ones(1,3);
normalsFaces = [eye(3) eye(3)]; % x y and z axes

points_line = zeros(3,6);
normals_line = zeros(3,6);
for i=1:6
    [pout, nout] = intersectionPlanePlane(normalsFaces(:,i),pointsFaces(:,i),slicing_plane.n, slicing_plane.p);
    points_line(:,i)=pout; % NaN????
    normals_line(:,i)=nout;
end

% find intersections between all lines
% get only points which are inside the bounding box
p_intersection=[];
npts = 0;
for i=1:5
    for j=(i+1):6
        Up =  normals_line(:,i);
        Uq =  normals_line(:,j);
        if (numel(find(isnan(Uq))) || numel(find(isnan(Up))) || Uq'*Up==1)
            % planes are somehow paralell
            continue;
        end
        P0 = points_line(:,i);
        Q0 = points_line(:,j);
        
        current_point =intersectionLineLine(Up,P0,Uq,Q0); % the closest point between both
        
        epsilon=10^(-5);
        if (current_point(1)>=(bounds(1)-epsilon) && current_point(1)<=(bounds(2)+epsilon) &&...
                current_point(2)>=(bounds(3)-epsilon) && current_point(2)<=(bounds(4)+epsilon) &&...
                current_point(3)>=(bounds(5)-epsilon) && current_point(3)<=(bounds(6)+epsilon))
            
            npts=npts+1;
            p_intersection(:,npts) = current_point ;
            intersecting_lines(npts,:)=[i j];
        end
    end
end


%% Calculate parameters of the 2D slice
if ~numel(p_intersection)
    % the slice is outside the volume
    disp('Error: the slice is outside the volume!');
    slice=NaN;
    return;
end
p_intersection_2D = M\ [p_intersection; ones(1,npts)];
% I think that the bounds_slice might be wrong...
bounds_slice = [min(p_intersection_2D,[],2); max(p_intersection_2D,[],2)];
bounds_slice = bounds_slice([1 5 2 6]);

if (spacing(1)==-1)
    spacing_slice = median(im.spacing)*ones(3,1);
else
    Mrot = M(1:3,1:3);
    spacing_slice  = Mrot\spacing;
    %spacing_slice = spacing_slice(1:numel(im.spacing));
end

size_slice = [ ceil((bounds_slice([2 4])-bounds_slice([1 3]))./spacing_slice(1:2)) ;1];
size_slice = max([size_slice ones(3,1)],[],2);
origin_slice= M*[ bounds_slice([1 3]) ; 0; 1]; % This is the origin in 3D
if thicknes>1
   size_slice(3) = 2*(thicknes-1)+1;
   origin_slice= M*[ bounds_slice([1 3]) ; 0-(thicknes-1)*spacing_slice(3); 1]; % This is the origin in 3D
else
    spacing_slice = [ spacing_slice(1:2) ; 1];   
end
origin_slice = origin_slice(1:3);
origin_slice2D=M\[origin_slice ; 1];
if (isa(im,'VectorImageType') || isfield(im,'datax'))
    slice = VectorImageType(size_slice , origin_slice, spacing_slice,M(1:3,1:3));
    slice2D = VectorImageType(size_slice([ 1 2]) , origin_slice2D([1 2]), spacing_slice([1 2]),eye(2));
else
    slice = ImageType(size_slice , origin_slice, spacing_slice,M(1:3,1:3));
    slice2D = ImageType(size_slice([1 2]) , origin_slice2D([1 2]), spacing_slice([1 2]),eye(2));
end

if (maskWithRoi)
    mask= ImageType(size_slice , [ bounds_slice([1 3]) ; 0], spacing_slice,eye(3));
    % convert roipoints to continuous indices
    roiindices = mask.GetContinuousIndex(roipoints);
    mask.data = roipoly( slice.data, roiindices(2,:), roiindices(1,:));
    
    
end

[x, y, z] = ndgrid(1:slice.size(1), 1:slice.size(2), 1:slice.size(3));
slicePositions = slice.GetPosition([x(:) y(:) z(:)]');
values = im.GetValue(slicePositions,interpolation);
values(values~=values)=0;
if (isa(im,'VectorImageType') || isfield(im,'datax'))
    if (maskWithRoi)
        slice.datax = reshape(values(1,:),slice.size').*mask.data;
        slice.datay = reshape(values(2,:),slice.size').*mask.data;
        slice.dataz = reshape(values(3,:),slice.size').*mask.data;
        slice.data = sqrt( slice.datax.^2 + slice.datay.^2 + slice.dataz.^2).*mask.data;
        
    else
        slice.datax = reshape(values(1,:),slice.size');
        slice.datay = reshape(values(2,:),slice.size');
        slice.dataz = reshape(values(3,:),slice.size');
        slice.data = sqrt( slice.datax.^2 + slice.datay.^2 + slice.dataz.^2);
    end
    
    % Convert to 2D
    vectors2D = M\[slice.datax(:) slice.datay(:) slice.dataz(:) 0*ones(numel(slice.datax),1)]';
    
    slice2D.datax(1:end) = vectors2D(1,:);
    slice2D.datax=slice2D.datax;
    slice2D.datay(1:end) = vectors2D(2,:);
    slice2D.datay = slice2D.datay;
    
    %slice2D.dataz(1:end) = dz;
    
else
    if (maskWithRoi)
        slice.data = reshape(values,slice.size').*mask.data;
    else
        slice.data = reshape(values,slice.size');
    end
    
    slice2D.data = slice.data;
    
end




%% Plots
if (doPlot)
    figure;
    hold on;
    s=200;
    lines_to_plot = unique(intersecting_lines(:));
    
    for i=lines_to_plot
        line([ points_line(1,i)-normals_line(1,i)*s ; points_line(1,i)+normals_line(1,i)*s] , ...
            [ points_line(2,i)-normals_line(2,i)*s ; points_line(2,i)+normals_line(2,i)*s] , ...
            [ points_line(3,i)-normals_line(3,i)*s ; points_line(3,i)+normals_line(3,i)*s] );
    end
    
    plot3(p_intersection(1,:),p_intersection(2,:),p_intersection(3,:),'*')
    
    % Draw bounds
    boundsPointsOrder = [1 2 4 3 1 5 6 2 6 8 4 8 7 3 7 5 ];
    line(pointsBounds(1,boundsPointsOrder ),pointsBounds(2,boundsPointsOrder ),pointsBounds(3,boundsPointsOrder ),'color',[1 0 0]);
    
    % draw plane
    
    m = planeMesh(slicing_plane.p, slicing_plane.n,'scale',s);
    trimesh(m.triangles,m.points(:,1),m.points(:,2),m.points(:,3) ,'faceColor',[0.4 0.4 0.4], 'edgeColor',[0.4 0.4 0.4]);
    
    hold off;
    axis equal;
    
end


end
