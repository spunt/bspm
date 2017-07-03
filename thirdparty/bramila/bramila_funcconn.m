function adj=bramila_funcconn(cfg)
    % add a check that the given connectivity measure exists and has all
    % that is needed to compute
    
    % add also the bit for the time of interests (toi)
	toi=1:size(cfg.roits,1);
	if(isfield(cfg,'toi') && ~isempty(cfg.toi))
		toi=cfg.toi;
		disp('Censoring time points')
	end

    switch cfg.conn_type
        case 'pearson'
            adj=corr(cfg.roits(toi,:),'type','pearson');
		case 'spearman'
            adj=corr(cfg.roits(toi,:),'type','spearman');
    end
    
