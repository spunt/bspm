function outm = surfaceFlatteningConstrained( inm_,varargin )
%SURFACEFLATTENING Flattens a surface mesh into a flat mesh
%   Using the method in:
%   M. Desbrun, M. Meyer, and P. Alliez, "Intrinsic Parametrization of
%   Surface Meshes", Eurographics 2002

lambda =[5 100];

do_natural_boundaries = false;

finalPos=[];
finalPosIdx=[];

for i=1:numel(varargin)
    if strcmp(varargin{i},'lambda')
        lambda = varargin{i+1};
    elseif strcmp(varargin{i},'finalPositions')
        finalPos = varargin{i+1};
        finalPosIdx = varargin{i+2};
    elseif strcmp(varargin{i},'naturalBoundary')
        do_natural_boundaries=true;
    end
end


% impose consistent orientation
%disp('Impose positive orientation. This will take a long while')
%inm.imposeConsistentOrientation();
%inm.imposeConsistentOrientationFlatMesh();
%disp('     orientation imposed.')
% find boundary
inm = MeshType(inm_);
id_border = inm.find_attribute('border');
if id_border>0
    b1 = inm.attributes(id_border).attribute_array;
else
    [b1,paths]=inm.detectBorders();
end


% Precalculate angles. Here we will calculate the angles alpha and delta as
% if the center of the ring was the first vertex of the triangle. The angle
% associated to that verted we call it mu (a23). Alpha = a12 and Delta =
% a13

DT_3D  = triangulation(inm.triangles,inm.points);
points = cat(3,DT_3D.Points(DT_3D.ConnectivityList(:,1),:),...
               DT_3D.Points(DT_3D.ConnectivityList(:,2),:),...
               DT_3D.Points(DT_3D.ConnectivityList(:,3),:));

sides12 = points(:,:,2)-points(:,:,1);
sides13 = points(:,:,3)-points(:,:,1);
sides23 = points(:,:,3)-points(:,:,2);

norms12=sum(sides12.^2,2); % length of each side squared
norms13=sum(sides13.^2,2); % length of each side squared
norms23=sum(sides23.^2,2); % length of each side squared

length12 = [norms12 norms23 norms13];
length12_extended = [length12 length12(:,1:2)];
angles = zeros(inm.ntriangles,3);

angles(:,1) = acosd((norms13 + norms23 - norms12)./(2*sqrt(norms13.*norms23))); % a12
angles(:,2) = acosd((norms12 + norms23 - norms13)./(2*sqrt(norms12.*norms23))); % a13
%a23 = acosd((norms12 + norms13 - norms23)./(2*sqrt(norms12.*norms13)));
angles(:,3) = 180-angles(:,1)-angles(:,2); % a23

angles_extended = [angles angles(:,[1 2])]; 

% beta = delta([2:end 1]);
% gamma= alpha([2:end 1]);


% points = cat(3,inm.points(inm.triangles(:,1),:),...
%     inm.points(inm.triangles(:,2),:),...
%     inm.points(inm.triangles(:,3),:));
% 
% sides12 = points(:,:,2)-points(:,:,1);
% sides13 = points(:,:,3)-points(:,:,1);
% sides23 = points(:,:,3)-points(:,:,2);
% 
% norms12=sum(sides12.^2,2); % ength of each side
% norms13=sum(sides13.^2,2); % ength of each side
% norms23=sum(sides23.^2,2); % ength of each side
% 
% alpha = acosd((norms13 + norms23 - norms12)./(2*sqrt(norms13.*norms23)));
% delta = acosd((norms12 + norms23 - norms13)./(2*sqrt(norms12.*norms23)));
% %a23 = acosd((norms12.^2 + norms13.^2 - norms23.^2)./(2*norms12.*norms13));



% Calculate the constant values intrinsic to the mesh

rows=[];
cols=[];

valsMA=[];
valsMX=[];

%indices_of_borders = find(b1);
indices_of_borders = find(b1==1);
not_borders = setdiff(1:inm.npoints,indices_of_borders);

if ( numel(finalPos))
    indices_of_borders = setdiff(indices_of_borders,finalPosIdx);
    not_borders = setdiff(not_borders ,finalPosIdx);
end




% get rings
rings = (DT_3D.vertexAttachments(not_borders(:)));


for i =1:numel(not_borders)
    
    % calculate the 1-ring neighbourhood, sorted
    
%      [ring , triangles1] = inm.GetRingOrdered(not_borders(i));
%      new_cols = [ triangles1(:,2)' not_borders(i)];
    
    

     ring = rings{i};
     triangles = DT_3D.ConnectivityList(ring,:);
     [rr,cc]=find(triangles==not_borders(i));
     triangles =[triangles triangles];
     [rr,idx2]=sort(rr);
     idx1 = sub2ind(size(triangles), rr,cc(idx2)+1);
     new_cols = [ triangles(idx1)' not_borders(i)];
     
     % Now we have to retrieve the angles and the length for the current
     % ring central vertex
     
     idx_alpha = sub2ind(size(angles_extended), ring(:),cc(idx2));
     idx_delta = sub2ind(size(angles_extended), ring(:),cc(idx2)+2);
     
     
     
    alpha = angles_extended(idx_alpha);
    delta = angles_extended(idx_delta);
    
    gamma= alpha([2:end 1]);  
    beta = delta([2:end 1]);

    
    length12sqared = length12_extended(idx_alpha);

    % Build matrices
    cols = [ cols new_cols];
    rows = [rows not_borders(i)*ones(1, numel(ring)+1)];
    inring_A = cotd(alpha)+cotd(beta);
    central_A = -sum(inring_A);
    
    valsMA = [valsMA inring_A' central_A];
    
    inring_X = (cotd(gamma)+cotd(delta))./length12sqared.^2; % This probably should not be squared, but squared works much better!
    central_X = -sum(inring_X);
    
    valsMX = [valsMX inring_X' central_X];
    
end


% for i =not_borders
%     
%     % calculate the 1-ring neighbourhood, sorted
%     [ring , triangles] = inm.GetRingOrdered(i);
%    
%     % Build matrices
%     rows = [rows i*ones(1, numel(ring)+1)];
%     cols = [ cols triangles(:,2)' i];
%     
%     inring_A = cotd(alpha(ring))+cotd(beta(ring));
%     central_A = -sum(inring_A);
%     
%     valsMA = [valsMA inring_A' central_A];
%     
%     inring_X = (cotd(gamma(ring))+cotd(delta(ring)))./norms12(ring).^2;
%     central_X = -sum(inring_X);
%     
%     valsMX = [valsMX inring_X' central_X];
%     
% end

if do_natural_boundaries
%     
%     rowsb=[];
%     colsb=[];
%     
%     valsMAb=[];
%     valsMXb=[];
%     
%     for i =indices_of_borders(:)'
%         
%         % calculate the 1-ring neighbourhood, sorted
%         %[~ , triangles] = inm.GetRingOrdered(not_borders(i));
%         [~ , triangles] = inm.GetRingOrdered(i);
%         
%         % calculate angles
%         
%         points = cat(3,inm.points(triangles(:,1),:),...
%             inm.points(triangles(:,2),:),...
%             inm.points(triangles(:,3),:));
%         
%         sides12 = points(:,:,2)-points(:,:,1);
%         sides13 = points(:,:,3)-points(:,:,1);
%         sides23 = points(:,:,3)-points(:,:,2);
%         norms12=sqrt(sum(sides12.^2,2)); % ength of each side
%         norms13=sqrt(sum(sides13.^2,2)); % ength of each side
%         norms23=sqrt(sum(sides23.^2,2)); % ength of each side
%         
%         alpha = acosd((norms13.^2 + norms23.^2 - norms12.^2)./(2*norms13.*norms23));
%         delta = acosd((norms12.^2 + norms23.^2 - norms13.^2)./(2*norms12.*norms23));
%         %a23 = acosd((norms12.^2 + norms13.^2 - norms23.^2)./(2*norms12.*norms13));
%         beta = delta([2:end 1]);
%         gamma= alpha([2:end 1]);
%         
%         % Build matrices
%         
%         
%         rowsb = [rowsb i*ones(1,size(triangles,1)+1)];
%         colsb = [ colsb triangles(:,2)' i];
%         
%         inring_A = cotd(alpha)+cotd(beta);
%         central_A = -sum(inring_A);
%         
%         valsMAb = [valsMAb inring_A' central_A];
%         
%         inring_X = (cotd(gamma)+cotd(delta))./norms12.^2;
%         central_X = -sum(inring_X);
%         
%         valsMXb = [valsMXb inring_X' central_X];
%         
%     end
%     
%     
%     rows = [rows  rowsb];
%     cols = [cols colsb];
%     
%     
%     MA = sparse(rows,cols,[valsMA valsMAb]*lambda(1),inm.npoints,inm.npoints);
%     MX = sparse(rows,cols,[valsMX valsMXb]*lambda(2),inm.npoints,inm.npoints);
%     
%     M = MA +MX;
%     
%     C = zeros(size(M,1),2);
%     d = zeros(numel(indices_of_borders));
%     for i=1:numel(indices_of_borders)
%         dist =inm.points(indices_of_borders,[1 2])-ones(numel(indices_of_borders),1)*inm.points(indices_of_borders(i),[1 2]);
%         d(:,i)=sum(dist.^2,2);
%     end
%     
%     [~,idx1D] = max(d(:));
%     [i1,i2]=ind2sub(size(d),idx1D);
%     
%     idx=[i1,i2];
%     %C(indices_of_borders(idx),:)=inm.points(indices_of_borders(idx),[1 2]);
%     C(indices_of_borders(idx),:)=inm.points(indices_of_borders(idx),[1 2]);
    
else
    
    % deal with borders
    rows = [rows  indices_of_borders(:)'];
    cols = [cols indices_of_borders(:)'];
    
    % deal with imposed positions
    
    rows = [rows  finalPosIdx(:)'];
    cols = [cols finalPosIdx(:)'];
    
    valsMA = [valsMA*lambda(1) zeros(1,numel(indices_of_borders)) zeros(1,numel(finalPosIdx))];
    valsMX = [valsMX*lambda(2) ones(1,numel(indices_of_borders)) ones(1,numel(finalPosIdx))];
    
    
    
    MA = sparse(rows,cols,valsMA,inm.npoints,inm.npoints);
    MX = sparse(rows,cols,valsMX,inm.npoints,inm.npoints);
    
    M = MA +MX;
    
    C = zeros(size(M,1),2);
    C(indices_of_borders,:)=inm.points(indices_of_borders,[1 2]);
    C(finalPosIdx,:)=finalPos(:,[1 2]);
    
end

% solve system
U = M\C;
outm = MeshType(inm);
outm.points = zeros(size(outm.points));
outm.points(:,[1 2])=U;

end

