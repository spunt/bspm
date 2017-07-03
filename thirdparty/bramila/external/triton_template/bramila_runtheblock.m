function bramila_runtheblock(id)
	% BRAMILA_RUNTHEBLOCK
	% This is a demo function for a potential parallel block. This function for example generates a random number
	% and stores it in a file uniquely identified by the block ID.
	% 
	% Add here your code for the specific block
	rng(id);
	a=rand;
	thispath=pwd;
	save([num2str(id) '.mat'],'a','thispath');
