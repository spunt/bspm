function [mea med perce sta]= bramila_scatterBars(data)
% this data is a matrix so that each column is a different group or condition and each row is a subject
% TODO: add colors, add options to plot/not-plot mean and/or median, add labels for the groups, checking for nans, estimate the distance with neighbours so that if distance is big enough then it can be plot in the centre

NG=size(data,2);
figure
for g=1:NG
	hold on
	thisdata=data(:,g);
	[AA BB]=sort(thisdata);
	shifts=repmat([-1 0 1]',length(BB),1);
	shifts=shifts(1:length(BB));
	thisdata=AA;
	mea=mean(thisdata);
	med=median(thisdata);	% I prefer median if there's enough data
	perce=prctile(thisdata,1:100);
	%thisdata=std(thisdata)./(thisdata+eps); % alternative, not in use
	plot(g+shifts/10,thisdata,'ko','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
	hold on
	plot([g-.15 g+.15],[mea mea],'r-','LineWidth',3)
   
	plot([g-.15 g+.15],[perce(75) perce(75)],'-','LineWidth',3,'color',[.5 .5 .5])
	plot([g-.15 g+.15],[perce(25) perce(25)],'-','LineWidth',3,'color',[.5 .5 .5])
	%plot(g,std(thisdata)/mea,'bo') % alternative for the coeff of variation
end

