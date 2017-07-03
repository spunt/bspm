function refdata = bramila_biopacFilter(refdata,filename)
% bramila_biopacFilter
% removes high apmplitudes and filter signas with FIR bandpass between 
% 0.05 Hz-1 Hz for breath 
% 0.3 Hz - 2 Hz for ppg
% saves figure to the *.fig file

breath=refdata{1}.data;
ppg=refdata{2}.data;
rawb=breath;
rawppg=ppg;
dt=refdata{1}.dt;

meanbreath=mean(breath);
stdbreath=std(breath);
f=breath>(meanbreath+3*stdbreath); 
breath(f)=(meanbreath+3*stdbreath+breath(f).*f(f)*0.05);
f=breath<(meanbreath-3*stdbreath);
breath(f)=(meanbreath-3*stdbreath+breath(f).*f(f)*0.05);

meanppg=mean(ppg);
stdppg=std(ppg);
f=ppg>(meanppg+3*stdppg); 
ppg(f)=(meanppg+3*stdppg+ppg(f).*f(f)*0.05);
f=ppg<(meanppg-3*stdppg);
ppg(f)=(meanppg-3*stdppg+ppg(f).*f(f)*0.05);

fs=1/dt;
nf=fs/2;

f1=0.05/nf;
f2=1/nf;
a=[1,0.5];
b=fir1(4000, [f1, f2]);
x=filter(b,a, [breath; zeros(2000,1)]);

f1ppg=0.3/nf;
f2ppg=2/nf;
a=[1,1];
b=fir1(2000, [f1ppg, f2ppg]);
y=filter(b,a, [ppg; zeros(1000,1)]);

refdata{1}.data=x(2001:end);
refdata{2}.data=y(1001:end);

h=figure('Visible','Off');
subplot(211)
plot([rawb, breath x(1:end-2000) refdata{1}.data])
title('breath')
legend('raw', 'amplitude normalized', 'filtered', 'filtered aligned = output')
subplot(212)
plot([rawppg, ppg y(1:end-1000) refdata{2}.data])
title('ppg')
legend('raw', 'amplitude normalized', 'filtered', 'filtered aligned = output')
saveas(h,[filename '_filtration.png'],'png');