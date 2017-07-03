function bramila_savepajek(cfg)

%	cfg.rois
%	cfg.adj
%	cfg.filename

%% it should run some input checks e.g. size of rois matches adj. Adj should be simmetric no neg values. Filename valid...

%% outputs the list of vertices
R = length(cfg.rois);
out = ['*Vertices ' num2str(R) '\n'];

for n = 1:R
    out=[out num2str(n) ' "' cfg.rois(n).label '" \n'];
end

% outputs adj matrix values that are > 0
out2='';
counter=0;
for r=1:R
	for c=(r+1):R
		val=cfg.adj(r,c);
		if(val>0);
			out2=[out2 num2str(r) ' ' num2str(c) ' ' num2str(val) '\n'];
			counter=counter+1;
		end
	end
end

out=[out '*Edges ' num2str(counter) '\n' out2];

fid=fopen(cfg.filename,'w');
fprintf(fid,out);
fclose(fid)
