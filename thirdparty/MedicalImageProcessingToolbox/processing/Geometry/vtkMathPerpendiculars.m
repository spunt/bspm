function [ y_, z_] = vtkMathPerpendiculars(x_,  theta)
% Given a unit vector 'x', find two other unit vectors 'y' and 'z' which
% which form an orthonormal set.
% Copied from vtkMath
% Modified to apply to any number of vectors
% x must be a 3 x N matrix

% do not compute vectors which are all 0

notAllZeros = find(x_(1,:)~=0 | x_(2,:)~=0 | x_(3,:)~=0);

x = x_(:,notAllZeros);

nvectors = size(x,2);

x2=x.^2;

r = sqrt(sum(x2,1));

% transpose the vector to avoid divide-by-zero error

dx = zeros(1,nvectors);
dy = zeros(1,nvectors);
dz = zeros(1,nvectors);

%  a = zeros(1,nvectors);
%  b = zeros(1,nvectors);
%  c = zeros(1,nvectors);

index1 = find(x2(1,:) > x2(2,:) & x2(1,:) > x2(3,:));
dx(index1)=1;
dy(index1)=2;
dz(index1)=3;
% a(index1)=x(1,:)./r;
% b(index1)=x(2,:)./r;
% c(index1)=x(3,:)./r;

if numel(index1)<nvectors
    index2 = find(x2(2,:) > x2(3,:));
        dx(index2)=2;
        dy(index2)=3;
        dz(index2)=1;
        
    if (numel(index1)+numel(index2))<nvectors
        index3 = setdiff(1:nvectors,union(index1,index2));
         dx(index3)=3;
        dy(index3)=1;
        dz(index3)=2;
%         a(index3)=x(3,:)./r;
%         b(index3)=x(1,:)./r;
%         c(index3)=x(2,:)./r;
    end
    
end

% convert dx, dy and dz to 1D indices

dx_ = sub2ind(size(x),dx,1:nvectors);
dy_ = sub2ind(size(x),dy,1:nvectors);
dz_ = sub2ind(size(x),dz,1:nvectors);


a=x(dx_)./r;    
b=x(dy_)./r;    
c=x(dz_)./r;    


% if (x2(1) > x2(2) && x2(1) > x2(3))
%     dx = 1; dy = 2; dz = 3;
% elseif (x2(2) > x2(3))
%     dx = 2; dy = 3; dz = 1;
% else
%     dx = 3; dy = 1; dz = 2;
% end
% 
% a = x(dx)/r;
% b = x(dy)/r;
% c = x(dz)/r;

tmp = sqrt(a.^2+c.^2);

y = zeros(3,nvectors);
z = zeros(3,nvectors);

if (theta ~= 0)
    sintheta = sin(theta);
    costheta = cos(theta);
    
    
    y(dx_) = (c*costheta - a.*b*sintheta)/tmp;
    y(dy_) = sintheta*tmp;
    y(dz_) = (-a*costheta - b.*c*sintheta)/tmp;
    
    
    
    z(dx_) = (-c*sintheta - a.*b*costheta)./tmp;
    z(dy_) = costheta*tmp;
    z(dz_) = (a*sintheta - b.*c*costheta)./tmp;
    
    
else
    
    y(dx_) = c./tmp;
    y(dy_) = 0;
    y(dz_) = -a./tmp;
    
    
    
    z(dx_) = -a.*b./tmp;
    z(dy_) = tmp;
    z(dz_) = -b.*c./tmp;
    
end


y_ = zeros(size(x_));
y_(:,notAllZeros)=y;
z_ = zeros(size(x_));
z_(:,notAllZeros)=z;

   
end