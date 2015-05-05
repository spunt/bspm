function h = plotpoints(p_,varargin)
% plots points using plot3
% p is a N x 3 array

    if isa(p_,'VectorImageType') || isa(p_,'ImageType')  || isfield(p_,'data')
        % count the number of images
        if numel(p_)>1
            % This is a set of images. plot the bounding box of each one
            % and its points
            mycolors=jet(numel(p_));
            hold on;
            for i=1:numel(p_)
               b=p_(i) .GetBounds();
               h(i,1)=plotboundingbox(b,varargin{1:end},'color',mycolors(i,:));
               h(i,2)=plotpoints(p_(i),'*','color',mycolors(i,:));
            end
            hold off;
            return;
        elseif numel(p_)==1
            % This is a single image plot the voxel centres
            p = p_.GetPosition(1:prod(p_.size))';
        end
    else
        p=p_;
    end
    
    if (nargin-1)
        h=plot3(p(:,1),p(:,2),p(:,3),varargin{1:end});
    else
        h=plot3(p(:,1),p(:,2),p(:,3));
    end
end