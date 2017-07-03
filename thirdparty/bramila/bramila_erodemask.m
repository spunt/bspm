function mask = bramila_erodemask(mask,erodemask)

% create voxel-free zone between two masks. Here the layer thickness is 1
% voxel, also counting diagonal neighbors.
% INPUT:
% mask = mask to be eroded (shrinks)
% erodemask = mask used for eroding operation (not modified)
% OUTPUT;
% mask = eroded mask (unchanged if two masks are already far separated)

MAX_NEAREST_NEIGHBORS = 0;

siz=size(mask);
x_max=siz(1);
y_max=siz(2);
z_max=siz(3);

[xx,yy,zz]=ind2sub(size(mask),find(mask));
N=length(xx);
for i=1:N
   x=(xx(i)-1):(xx(i)+1);
   y=(yy(i)-1):(yy(i)+1);
   z=(zz(i)-1):(zz(i)+1);
   
   x(x<1 | x>x_max)=[];
   y(y<1 | y>y_max)=[];
   z(z<1 | z>z_max)=[];
   
   val = sum(sum(sum(erodemask(x,y,z))));
   if val > MAX_NEAREST_NEIGHBORS % remove if too many voxels inside cube
       mask(xx(i),yy(i),zz(i))=0;
   end
end

% CONN implementation:
% X0=in_mask;
% idx1=find(X0(:));
% [idxx,idxy,idxz]=ind2sub(size(X0),idx1);
% idxt=find(idxx>ERODE&idxx<size(X0,1)+1-ERODE&idxy>ERODE&idxy<size(X0,2)+1-ERODE&idxz>ERODE&idxz<size(X0,3)+1-ERODE);
% for n1=1:length(idxt),
%     if (sum(sum(sum(X0(idxx(idxt(n1))+(-ERODE:ERODE),idxy(idxt(n1))+(-ERODE:ERODE),idxz(idxt(n1))+(-ERODE:ERODE))<.5,3),2),1))>1,
%         idxt(n1)=0;
%     end;
% end
% idxt=idxt(idxt>0);
% idx1=idx1(idxt);
% out_mask=zeros(size(X0));
% out_mask(idx1)=1;

end

