function reg_cfg = bramila_makemotionregress(cfg)

    X_motion = load(cfg.motionparam);
    
    labels = {'X','Y','Z','pitch','yaw','roll'};

	prepro_suite='fsl-fs'; 
    if(isfield(cfg,'prepro_suite'))
        prepro_suite = cfg.prepro_suite;
    end

    radius=50; % default radius
    if(isfield(cfg,'radius'))
        radius = cfg.radius;
    end

    if(strcmp(prepro_suite,'fsl-fs'))
		temp = X_motion;
		X_motion=[temp(:,4:6) temp(:,1:3)]; % in fsl we swap rotations and translations so that we match our labels on top
	end

    X_motion(:,4:6)=X_motion(:,4:6)*radius;	% they are already in radians

    X_motion_der = [];
    labels_der = [];
    for i=1:cfg.mot_derivatives
        X_motion_der=cat(2,X_motion_der,bramila_derivative(X_motion,i));       
        labels_der = cat(2,labels_der,get_deriv_labels(labels,i));
    end
    
    % detrend motion + derivatives before building regressors
    temp_cfg.TR=cfg.TR;
    temp_cfg.detrend_type=cfg.detrend_type; % must be the same as applied to BOLD data
    temp_cfg.write=0;
    
    temp_cfg.vol=X_motion;
    X_motion = bramila_detrend(temp_cfg);
    if ~isempty(X_motion_der)
        temp_cfg.vol=X_motion_der;
        X_motion_der = bramila_detrend(temp_cfg);
    end
    
    % build motion design matrix        
    if strcmp(cfg.motion_reg_type,'standard') % standard (LAB CLASSIC)
        % no additional regressors needed
    elseif strcmp(cfg.motion_reg_type,'friston') % Friston style with squares and shifts
        frist1=circshift(X_motion,[1,0]);
        frist1(1,:)=0;                
        X_motion=[X_motion,X_motion.^2,frist1,frist1.^2];
        labels = {'X','Y','Z','pitch','yaw','roll',...
            'X^2','Y^2','Z^2','pitch^2','yaw^2','roll^2',...
            'X (+1TR)','Y (+1TR)','Z (+1TR)','pitch (+1TR)','yaw (+1TR)','roll (+1TR)',...
            'X^2 (+1TR)','Y^2 (+1TR)','Z^2 (+1TR)','pitch^2 (+1TR)','yaw^2 (+1TR)','roll^2 (+1TR)'};
    elseif strcmp(cfg.motion_reg_type,'volterra') % total 36 regressors
        frist1=circshift(X_motion,[1,0]);
        frist1(1,:)=0;
        frist2=circshift(X_motion,[2,0]);
        frist2(1:2,:)=0;        
        X_motion=[X_motion,X_motion.^2,frist1,frist1.^2,frist2,frist2.^2];
        labels = {'X','Y','Z','pitch','yaw','roll',...
            'X^2','Y^2','Z^2','pitch^2','yaw^2','roll^2',...
            'X (+1TR)','Y (+1TR)','Z (+1TR)','pitch (+1TR)','yaw (+1TR)','roll (+1TR)',...
            'X^2 (+1TR)','Y^2 (+1TR)','Z^2 (+1TR)','pitch^2 (+1TR)','yaw^2 (+1TR)','roll^2 (+1TR)',...
            'X (+2TR)','Y (+2TR)','Z (+2TR)','pitch (+2TR)','yaw (+2TR)','roll (+2TR)',...
            'X^2 (+2TR)','Y^2 (+2TR)','Z^2 (+2TR)','pitch^2 (+2TR)','yaw^2 (+2TR)','roll^2 (+2TR)'};        
    else
        error('Unknown motion regression type!')
    end

    X_total = [X_motion,X_motion_der]; % add derivatives (if present)
    labels_total = cat(2,labels,labels_der);

    reg_cfg.reg = X_total;
    reg_cfg.vol = cfg.vol;
    reg_cfg.labels = labels_total;

end

function labels = get_deriv_labels(labels,order)

    dot = repmat('''',1,order);
    for i=1:length(labels)
        labels{i}=[labels{i},dot];
    end
end

% POWERS et al.:
%         switch switches.motionestimates
%             case 0 %
%                 QC(i).mvmregs=[];
%                 QC(i).mvmlabels={''};
%             case 1 % R,R`                   LAB CLASSIC
%                 QC(i).mvmregs=[QC(i).DTMVM QC(i).ddtDTMVM];
%                 QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','X`','Y`','Z`','pitch`','yaw`','roll`'};
%             case 2 % R,R^2,R-1,R-1^2       FRISTON
%                 frist1=circshift(QC(i).DTMVM,[1 0]);
%                 frist1(1,:)=0;
%                 QC(i).mvmregs=[QC(i).DTMVM (QC(i).DTMVM.^2) frist1 frist1.^2 ];
%                 QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','sqrX','sqrY','sqrZ','sqrpitch','sqryaw','sqrroll`','Xt-1','Yt-1','Zt-1','pitcht-1','yawt-1','rollt-1`','sqrXt-1','sqrYt-1','sqrZt-1','sqrpitcht-1','sqryawt-1','sqrrollt-1`'};
%             case 20 % R,R`,12               CONTROL for 24 parameters
%                 QC(i).mvmregs=[QC(i).DTMVM QC(i).ddtDTMVM rand(size(QC(i).DTMVM,1),12) ];
%                 QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','X`','Y`','Z`','pitch`','yaw`','roll`','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12'};
%             case 3
%                 frist=[QC(i).DTMVM];
%                 frist1=circshift(QC(i).DTMVM,[1 0]);
%                 frist1(1,:)=0;
%                 frist2=circshift(QC(i).DTMVM,[2 0]);
%                 frist2(1:2,:)=0;
%                 QC(i).mvmregs=[QC(i).DTMVM (QC(i).DTMVM.^2) frist1 frist1.^2 frist2 frist2.^2 ];
%                 QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','sqrX','sqrY','sqrZ','sqrpitch','sqryaw','sqrroll`','Xt-1','Yt-1','Zt-1','pitcht-1','yawt-1','rollt-1`','sqrXt-1','sqrYt-1','sqrZt-1','sqrpitcht-1','sqryawt-1','sqrrollt-1`','Xt-2','Yt-2','Zt-2','pitcht-2','yawt-2','rollt-2`','sqrXt-2','sqrYt-2','sqrZt-2','sqrpitcht-2','sqryawt-2','sqrrollt-2`'};
%         end
