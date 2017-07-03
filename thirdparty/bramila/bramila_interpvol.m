function bramila_interpvol(cfg);
%
% BRAMILA_INTERPVOL - Interpolates a volume e.g. from 6mm to 2mm space, useful to visualize our results with 6mm nodes
%   - Usage:
%   fwd = bramila_interpvol(cfg)
%   - Input:
%   cfg is a struct with following parameters
%		cfg.rois = rois definition;
%		cfg.blacklist = list of nodes to exclude
%		cfg.vals = values to plot for each node (column vector of length == length(rois)). 
%			Note: it can contain more than one column, output will be saved as 4-D image
%		cfg.datasize = the standard template output ([91 109 91] in our cases)
%		cfg.box - 0 (put > 0 if you want a smoothing filter)
%		cfg.interptype = linear, cubic, spline or none
%		cfg.filename = a filename.nii where the data will be saved				
%   - Output:
%       No output, only writes to file. Note that the file has originator issues.
%   - Notes:
%		Add support for bramila_savevolume to fix originator issues. Add better user input testing.
%	Last edit: EG 2014-09-10

    rois=cfg.rois;
    blacklist=cfg.blacklist;
    R=length(rois);
    vals=cfg.vals;
    vals(blacklist,:)=0;

    if(R ~= size(vals,1))
        error('mismatch in the length')
    end
	
	% more than one image at once
	Nvals=size(vals,2);
	allvols=zeros([cfg.datasize(1:3) Nvals]);
	for currval=1:Nvals
            mask=zeros(cfg.datasize(1:3));
		vol=zeros(cfg.datasize(1:3));
		%--> load(psess.groupmask) add this at the end!!!
        
		%vol(find(groupmask==1))=NaN;
		%compare=vol;
		X=zeros(R,1);
		Y=X;
		Z=X;
		Xi=[];
		Yi=[];
		Zi=[];
		outnone=zeros(cfg.datasize(1:3));
		for r=1:R
			c=rois(r).centroid;
			X(r)=c(1);
			Y(r)=c(2);
			Z(r)=c(3);
			
			
		    map=rois(r).map;
            for m=1:size(map,1)
                    x=map(m,1);
		        y=map(m,2);
		        z=map(m,3);

		        mask(x,y,z)=1;
				outnone(x,y,z)=vals(r,currval);
		    end
			
			vol(c(1),c(2),c(3))=vals(r,currval);
            
		end
		
		if(strcmp(cfg.interptype,'none')==1)
			voli=outnone;
		else
			[gX,gY,gZ] = meshgrid(2:3:109,2:3:91,2:3:91);
			[bgX, bgY,bgZ] = meshgrid(1:cfg.datasize(2),1:cfg.datasize(1),1:cfg.datasize(3));
			voli=interp3(gX,gY,gZ,vol(2:3:91,2:3:109,2:3:91),bgX, bgY,bgZ,cfg.interptype);
		end
		voli(find(isnan(voli)))=0;
        amp=1;
        if(cfg.box>0)
            voli=smooth3(voli,'gaussian',[ 3 3 3]);     % gaussian leaves the original data a bit more intact...
            %voli=smooth3(voli,'box',cfg.box);
            % we need to scale the smoothed data
            ids=find(vol>max([median(vals(:,currval)) 0]));
            amp=mean(voli(ids)./vol(ids));
            
            if(isnan(amp))
                save workspace_debug
                error('its nan!')
            end
        end
        
        
		%allvols(:,:,:,currval)=mask.*voli/amp;
        allvols(:,:,:,currval)=voli/amp;
        
	end
	save_nii(make_nii(allvols),cfg.filename);
    
    
    %save_nii(make_nii(vol),cfg.filename);
	%voli=smooth3(voli,'gaussian',[ 3 3 3]);
    %save_nii(make_nii(voli),'testinterp_linear_smoo.nii');
    
%    voli=interp3(gX,gY,gZ,vol(2:3:91,2:3:109,2:3:91),bgX, bgY,bgZ,'spline');
%    save_nii(make_nii(voli),'testinterp_spline.nii');
%    voli=smooth3(voli,'gaussian',[ 3 3 3]);
%    save_nii(make_nii(voli),'testinterp_spline_smoo.nii');
       
    
    
    %voli=inpaint_nans3(vol);
    
    
    %save_nii(make_nii(compare),'compareinterp.nii');
    
