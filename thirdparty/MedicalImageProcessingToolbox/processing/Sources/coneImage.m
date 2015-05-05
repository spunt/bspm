% test_generateCone
function out = coneImage(points_axis, points_side, ref_im)


%% parameters

out = ImageType(ref_im);
 th = sqrt(sum(ref_im.spacing.^2))/2;
%th=48;
npieces=1;
%% get the apex of the cone

Up =  (points_axis(:,2)-points_axis(:,1))/norm(points_axis(:,2)-points_axis(:,1));
Uq =  (points_side(:,2)-points_side(:,1))/norm(points_side(:,2)-points_side(:,1));

P0 = points_axis(:,1);
Q0 = points_side(:,1);

b = -1*[(P0-Q0)'*Up; (P0-Q0)'*Uq];
A = [ Up'*Up -Uq'*Up;  Up'*Uq -Uq'*Uq];

lambda = A\b;

P1 = P0 + lambda(1)*Up;
%Q1 = Q0 + lambda(2)*Uq;

apex = P1;

%% Get distance to the cone

[ x y z] = ndgrid(1:ref_im.size(1), 1:ref_im.size(2), 1:ref_im.size(3));
positions = ref_im.GetPosition([x(:) y(:) z(:)]');

vectors = positions - apex*ones(1,size(positions,2));
vectors_norm = sqrt(vectors(1,:).^2+vectors(2,:).^2+vectors(3,:).^2);

vectors(1,:)=vectors(1,:)./vectors_norm;
vectors(2,:)=vectors(2,:)./vectors_norm;
vectors(3,:)=vectors(3,:)./vectors_norm;

angle_with_coneAxis = acos(abs(dot(vectors,(Up*ones(1,size(vectors,2))),1)));

cone_angle= acos(abs(dot(Uq,Up)));

distance_from_point_to_cone =  vectors_norm .* sin(cone_angle*ones(1,size(angle_with_coneAxis,2)) -angle_with_coneAxis);
%distance_from_point_to_cone =  abs(angle_with_coneAxis-cone_angle*ones(1,size(angle_with_coneAxis,2)));
%distance_from_point_to_cone(distance_from_point_to_cone<0)=0;

out.data = reshape(abs(distance_from_point_to_cone)<th,out.size');
%out.data = reshape(distance_from_point_to_cone,out.size');
%out.data = reshape(angle_with_coneAxis*180/pi,out.size');

% out2=VectorImageType(ref_im);
% out2.datax = reshape(vectors(1,:),out2.size');
% out2.datay = reshape(vectors(2,:),out2.size');
% out2.dataz = reshape(vectors(3,:),out2.size');




end
