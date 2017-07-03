function volout=bramila_tsnr(volin);

% Usage:
% 	tSNR = bramila_tsnr(data)
%		data = a four dimensional volume

M=mean(volin,4);
S=std(volin,0,4);
volout=M./S;
