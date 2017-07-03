function bramila_savematrix(cfg,mat)

if ~exist([cfg.inpath,cfg.separator,'bramila'],'dir')
    mkdir([cfg.inpath,cfg.separator,'bramila']);
end

mat=double(mat);
%if max(max(abs(mat-mat')))<1e-13
%   mat = sparse(triu(mat,1));
%end

adj=mat;
save([cfg.outpath,cfg.separator,'bramila',cfg.separator,cfg.fileID,'_adj.mat'],'adj','-v7.3');

end
