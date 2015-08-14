function ret = ipctb_spm_ov_rgb(varargin)
% RGB overlays
% A shorthand to overlaying the absolute value of three different images
% onto a displayed image in colours red, green and blue. The overlay images
% are optionally masked and multiplied with a scaling image. The displayed
% overlay images are the absolute value of the given overlays.
%
% This routine is a plugin to ipctb_spm_orthviews for SPM5. For general help about
% ipctb_spm_orthviews and plugins type
%             help ipctb_spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% @(#) $Id: ipctb_spm_ov_rgb.m 591 2006-08-14 11:06:49Z volkmar $

rev = '$Revision: 591 $';

global st;
if isempty(st)
  error('rgb: This routine can only be called as a plugin for ipctb_spm_orthviews!');
end;

if nargin < 2
  error('rgb: Wrong number of arguments. Usage: ipctb_spm_orthviews(''rgb'', cmd, volhandle, varargin)');
end;

cmd = lower(varargin{1});
volhandle = varargin{2};

switch cmd
  
  %-------------------------------------------------------------------------
  % Context menu and callbacks
  case 'context_menu'  
    item0 = uimenu(varargin{3}, 'Label', 'RGB overlays');
    item1 = uimenu(item0, 'Label', 'Add', 'Callback', ...
	['feval(''ipctb_spm_ov_rgb'',''context_init'', ', ...
	  num2str(volhandle), ');'], 'Tag', ['RGB_0_', num2str(volhandle)]);
    item1 = uimenu(item0, 'Label', 'Help', 'Callback', ...
	  ['feval(''ipctb_spm_help'',''' mfilename ''');']);
    
  case 'context_init'
    Finter = ipctb_spm_figure('FindWin', 'Interactive');
    ipctb_spm_input('!DeleteInputObj',Finter);
    [Vqfnames sts] = ipctb_spm_select(3, 'image',...
                                'Components of 1st eigenvector', ...
                                [], pwd, 'evec1.*');
    if ~sts return; end;
    Vq = ipctb_spm_vol(Vqfnames);
    [Vfafname sts] = ipctb_spm_select([0 1],'image','FA image (optional)', [], ...
                                pwd, 'fa.*'); 
    Vfa = ipctb_spm_vol(Vfafname);
    [Vmaskfname sts] = ipctb_spm_select([0 1],'image','Mask image (optional)');
    Vmask = ipctb_spm_vol(Vmaskfname);
    ipctb_spm('pointer','watch');
    Vamq = rmfield(Vq,'private');
    for k=1:3
      [p n e v]=fileparts(Vq(k).fname);
      sel = 2*isempty(Vmask)+isempty(Vfa);
      switch(sel)
	case 0, %both Vmask and Vfa set
	  Vamq(k).fname=fullfile(p,['abs_msk_fa_' n e v]);
	  ipctb_spm_imcalc([Vq(k) Vfa Vmask],Vamq(k),'abs(i1.*i2.*i3)',{[],1,[]});
	case 1, %only Vmask set
	  Vamq(k).fname=fullfile(p,['abs_msk_' n e v]);
	  ipctb_spm_imcalc([Vq(k) Vmask],Vamq(k),'abs(i1.*i2)',{[],1,[]});
	case 2, %only Vfa set
	  Vamq(k).fname=fullfile(p,['abs_fa_' n e v]);
	  ipctb_spm_imcalc([Vq(k) Vfa],Vamq(k),'abs(i1.*i2)',{[],1,[]});
	case 3, %nothing set
	  Vamq(k).fname=fullfile(p,['abs_' n e v]);
	  ipctb_spm_imcalc(Vq(k),Vamq(k),'abs(i1)',{[],1,[]});
      end;
    end;
    ipctb_spm_orthviews('addcolouredimage',volhandle,Vamq(1).fname,[1 0 0]);
    ipctb_spm_orthviews('addcolouredimage',volhandle,Vamq(2).fname,[0 1 0]);
    ipctb_spm_orthviews('addcolouredimage',volhandle,Vamq(3).fname,[0 0 1]);
    ipctb_spm_orthviews('redraw');
    
    ipctb_spm_input('!DeleteInputObj',Finter);
 case 'redraw'
  % Do nothing
  otherwise    
    fprintf('ipctb_spm_orthviews(''rgb'', ...): Unknown action %s', cmd);
end;

ipctb_spm('pointer','arrow');
