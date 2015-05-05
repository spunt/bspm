function h = plotAxis(direction, varargin)
% plots 3 vectors with xyz, (columns of the 3x3 matrix direction)
%

scale=1;
currentArg=1;
ls='-';
lw=2;
if numel(direction)==9
    origin = varargin{1};
    currentArg=2;
else
    % matrix is 4x4
    origin = direction(1:3,4);
end



for i=currentArg:size(varargin,2)
    if  (strcmp(varargin{i},'scale'))
        scale = varargin{i+1};
        i=i+1;
    elseif  (strcmp(varargin{i},'LineStyle'))
        ls= varargin{i+1};
        i=i+1;
     elseif  (strcmp(varargin{i},'LineWidth'))
        lw= varargin{i+1};
        i=i+1;
    end
    
end

hold on;
h.h1 = quiver3(origin(1),origin(2),origin(3),direction(1,1),direction(2,1),direction(3,1),scale,'Color',[1 0 0],'LineWidth',lw,'LineStyle',ls);
h.h2 = quiver3(origin(1),origin(2),origin(3),direction(1,2),direction(2,2),direction(3,2),scale,'Color',[0 1 0],'LineWidth',lw,'LineStyle',ls);
h.h3 = quiver3(origin(1),origin(2),origin(3),direction(1,3),direction(2,3),direction(3,3),scale,'Color',[0 0 1],'LineWidth',lw,'LineStyle',ls);




end