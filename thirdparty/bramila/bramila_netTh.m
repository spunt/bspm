function [idsout, MST]=bramila_netTh(A,th,MST)
%	A = adj matrix
%	th = vector of percetages of link density
R=size(A,1);
L=(R^2-R)/2;
ids=find(triu(ones(size(A)),1));
[vals edgt]=sort(A(ids),'descend');
edg=ids(edgt);

if(size(MST,1)>1)
	MST=sign(MST);
else
	%% build mst
	GCC=0;
	MST=zeros(size(A));
	counter=1;
	oldcosiz=-1;
	while(GCC<R)
		temp=triu(MST,1);
		temp(edg(counter))=1;
		temp=temp+temp';
		[co,cosiz]=get_components(temp);
		if(max(cosiz)>oldcosiz)
			MST=temp;
			oldcosiz=max(cosiz);
			GCC=max(cosiz);
		end
		counter=counter+1;
		fprintf([num2str(counter) '..'])
	end

end
idsMST=find(triu(MST,1));
size(idsMST)
edg=setdiff(edg,idsMST,'stable');
edg=[idsMST;edg];

tt=round(th*L/100);
for t=1:length(tt)
	idsout{t}=edg(1:tt(t));
	%temp=zeros(size(A));
	%temp(idsout{t})=1;
	%temp=temp+temp';
	%[co,cosiz]=get_components(temp);
	%if(max(cosiz)<R)
%		disp([num2str(th(t)) ' would need an MST']);
%	end
end


