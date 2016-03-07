function drawlines(x,color,range)
% Drawlines(x,color,range)
% range is a ymin-ymax vector (set to size of axis as default)
% draws x-seperation line into a graph at the places of an x-vector
if (nargin<3)
      range=get(gca,'YLim');
end;
[r,c]=size(x);
if(r>c)
   x=x';
end;
x1=[x;x];
y=[ones(1,length(x));ones(1,length(x))];
y(1,:)=y(1,:).*range(1);
y(2,:)=y(2,:).*range(2);
line(x1,y,'Color',color);
