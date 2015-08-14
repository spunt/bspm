function ret = ipctb_spm_ov_movie(varargin)
% Movie tool - plugin for ipctb_spm_orthviews
%
% This plugin allows an automatic "fly-through" through all displayed
% volumes. Apart from pre-defined trajectories along the x-, y- and z-axis,
% resp., it is possible to define custom start and end points (in mm) for
% oblique trajectories.
%
% This routine is a plugin to ipctb_spm_orthviews for SPM5. For general help about
% ipctb_spm_orthviews and plugins type
%             help ipctb_spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% @(#) $Id: ipctb_spm_ov_movie.m 591 2006-08-14 11:06:49Z volkmar $

global st;
if isempty(st)
  error('movie: This routine can only be called as a plugin for ipctb_spm_orthviews!');
end;

if nargin < 2
  error('movie: Wrong number of arguments. Usage: ipctb_spm_orthviews(''roi'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd
  
  %-------------------------------------------------------------------------
  % Context menu and callbacks
  case 'context_menu'  
    item0 = uimenu(varargin{3}, 'Label', 'Movie tool');
    item1 = uimenu(item0, 'Label', 'Run', 'Callback', ...
	['feval(''ipctb_spm_ov_movie'',''context_init'', ', ...
	  num2str(volhandle), ');'], 'Tag', ['MOVIE_0_', num2str(volhandle)]);
    item1 = uimenu(item0, 'Label', 'Help', 'Callback', ...
	  ['feval(''ipctb_spm_help'',''' mfilename ''');']);

  case 'context_init'
    Finter = ipctb_spm_figure('FindWin', 'Interactive');
    opos=ipctb_spm_orthviews('pos');
    ipctb_spm_input('!DeleteInputObj',Finter);
    dir=logical(cell2mat(ipctb_spm_input('Select movie direction', '!+1', 'b', 'x|y|z|custom', ...
	{[1 0 0], [0 1 0], [0 0 1], 0}, 1)));
    if all(dir==0)
      mstart=ipctb_spm_input('First point (mm)', '!+1', 'e', num2str(opos'), [3 1]);
      mend  =ipctb_spm_input('Final point (mm)', '!+1', 'e', num2str(opos'), [3 1]);
    else
      mstart=opos;
      mend=opos;
      bb = st.Space*[st.bb'; 1 1];
      dirs='XYZ';
      tmp=ipctb_spm_input([dirs(dir) ' intervall (mm)'], '!+1', 'e', ...
	  num2str(bb(dir,:), '%.1f %.1f'), 2);
      mstart(dir)=tmp(1);
      mend(dir)=tmp(2);
    end;
    ds=ipctb_spm_input('Step size (mm)', '!+1', 'e', '1', 1);
    d=mend-mstart;
    l=sqrt(d'*d);
    d=d./l;
    for k=0:ds:l
      ipctb_spm_orthviews('reposition', mstart+k*d);
    end;
    ipctb_spm_orthviews('reposition', opos);
    ipctb_spm_input('!DeleteInputObj',Finter);
  otherwise    
    fprintf('ipctb_spm_orthviews(''movie'', ...): Unknown action %s', cmd);
end;


ipctb_spm('pointer','arrow');
