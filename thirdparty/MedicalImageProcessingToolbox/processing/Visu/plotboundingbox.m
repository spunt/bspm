function h = plotboundingbox(bounds,varargin)
% bounds is a column vector


ndimensions = numel(bounds)/2;

m = zeros(ndimensions,2);

for i=1:ndimensions
    m(i,1)=bounds(2*(i-1)+1);
    m(i,2)=bounds(2*(i-1)+2);
end

% create the 2^ndimensions corners
str='';
str2='';
str3='';
for i=1:ndimensions
   str=[str '[1 2],']; 
   str2=[str2 'x' num2str(i) ','];
   str3=[str3 'x' num2str(i) '(:) '];
end
eval([ '[' str2(1:end-1) ']=ndgrid(' str(1:end-1) '); indices = [' str3(1:end-1) '];'  ]);

corners = zeros(2^ndimensions,1);

for r=1:size(indices,1)
    for i=1:ndimensions
        corners(r,i)=m(i,indices(r,i));
    end
end

    if ndimensions==3
        order = [1 3 4 2 1 5 7 8 6 5 6 2 6 8 4 8 7 3];
        if (nargin-1)
            h=plot3(corners(order,1),corners(order,2),corners(order,3),varargin{1:end});
        else
            h=plot3(corners(order,1),corners(order,2),corners(order,3));
        end
    elseif ndimensions==2
        order = [1 3 4 2 1];
        if (nargin-1)
            h=plot(corners(order,1),corners(order,2),varargin{1:end});
        else
            h=plot(corners(order,1),corners(order,2));
        end
    else
        disp('plotboundingbox only works for 2D or 3D')
    end
    
end