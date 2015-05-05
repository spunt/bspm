function h = plotpoints2(p_,varargin)
% plots points using plot3
% p is a N x 3 array

    if isa(p_,'VectorImageType') || isa(p_,'ImageType')  || isfield(p_,'data')
        % plot the voxel centres
        p = p_.GetPosition(1:prod(p_.size))';
    else
        p=p_;
    end
    
    if (nargin-1)
        h=plot(p(:,1),p(:,2),varargin{1:end});
    else
        h=plot(p(:,1),p(:,2));
    end
end