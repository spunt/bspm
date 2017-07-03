function [x,y,z] = bramila_MNI(cfg)

% BRAMILA_MNI - Converts 2mm MNI coordinates into matrix indeces or viceversa.
% - Usage:
% 	cfg.type='MNI' or 'space' - the type of the input
%	cfg.coordinates = a matrix with 3 column, one row per location
%	cfg.imgsize = [91 109 91] (fsl) or [79 95 68] (spm) 
% [x,y,z] = bramila_MNI2space(cfg);
%
% EG 2014-01-14
% Notes: 	add support for more millimeters sizes 
% 			see also: http://brainmap.wustl.edu/help/mapper.html

switch cfg.type
	case 'MNI'
		xMNI = cfg.coordinates(:,1);
		yMNI = cfg.coordinates(:,2);
		zMNI = cfg.coordinates(:,3);

		if (cfg.imgsize == [91 109 91])
			x = xMNI/2 + 46;
			y = 108*(yMNI + 126)/216 +1;
			z=90*(zMNI + 72)/180 +1;
		end
		if(cfg.imgsize == [79 95 68])
			x = (xMNI+78)/2+1;
			y = 94*(yMNI + 112)/(76+112)+1;
			z=67*(zMNI+50)/(84+50)+1;
		end
	case 'space'
		xsp = cfg.coordinates(:,1);
        ysp = cfg.coordinates(:,2);
        zsp = cfg.coordinates(:,3);

		if (cfg.imgsize == [91 109 91])
			x = 2*((xsp-1)-45);
		    y = 216*(ysp-1)/108 - 126;
    		z = 180*(zsp-1)/90 - 72;
		end
		if (cfg.imgsize == [79 95 68])
			x = 2*(xsp-1)-78;
			y = (ysp -1)*(76+112)/94 - 112;
			z = (zsp -1)*(84+50)/67 - 50;
		end

	otherwise 
		disp(['Type uknown ' cfg.type])
end
