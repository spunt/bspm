function bramila_plotConn(tnet,rois,titlelabel,imagesclimits,Nwhites);
	tempmap=flipud(cbrewer('div','RdBu',11));

	map=[tempmap(1:5,:); ones(Nwhites,3); tempmap(7:end,:)];
	%map(6,:)=[1 1 1];
	

	R=length(rois);
	for r=1:R
		labels{r}=rois(r).label;
	end

	imagesc(tril(tnet),imagesclimits);
	% store the final tnet eventually
	hold on
	% plot the diagonal
	colormap(map)

	colorbar('SouthOutside')

	axis square

	for cc=1:R
		plot([cc-.5 cc-.5],[cc-.5 cc+.5],'Color',[.5 .5 .5]);
		plot([cc-.5 cc+.5],[cc+.5 cc+.5],'Color',[.5 .5 .5]);
		plot([cc+.5 cc+.5],[cc+.5 cc-.5],'Color',[.5 .5 .5]);
		plot([cc+.5 cc-.5],[cc-.5 cc-.5],'Color',[.5 .5 .5]);

		% full lines
		plot([0 cc]+.5,[cc cc]+.5,'Color',[.5 .5 .5]);
		plot([cc cc]-.5,[cc R]+.5,'Color',[.5 .5 .5]);


		t=text(cc+.5,cc-.5,labels{cc},'FontSize',7,'Interpreter','none','Rotation',45,'FontName','Arial');
		set(gca,'YTick',[1:119]);
		set(gca,'YTickLabel',labels);
		%set(gca,'FontSize',6)
		set(gcf,'Color',[1 1 1])
		axis off
		box off
		title(titlelabel);
	end

