function varargout=osdisp(th0,thhats,nl,avhs,Fisher,truecov)
% OSDISP(th0,thhats,nl,avhs,Fisher,truecov)
% OSDISP(th0,params)
% [str0,str1,str2,str3]=OSDISP(...)
%
% Displays some messages to the command screen for MLEOS, MLEROS, MLEROS0
% and also SIMULOS, SIMULROS, SIMULROS0
%
% INPUT:
%
% th0        True parameter vector
% thhats     Estimated parameter vector, OR
% params     A structure with the fixed parameter settings
% nl         Number of experiments over which the average Hessian is reported
% avhs       Average Hessian matrix
% Fisher     Fisher matrix
% truecov    The theoretical parameter covariance
%
% OUTPUT:
%
% The strings used 
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2014

% The necessary strings for formatting
str0='%27s';
str1='%12.5g ';
str2='%12i ';
str3='%12s ';

% Replicate to the size of the parameter vector
str1s=repmat(str1,size(th0));
str2s=repmat(str2,size(th0));

if length(th0)==6
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str3,size(th0))),...
	       'True theta','D','f2','r','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str1s),' ',th0))
elseif length(th0)==5
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str3,size(th0))),...
	       'True theta','D','f2','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str1s),' ',th0))
elseif length(th0)==3
  disp(sprintf(sprintf('\n%s : %s ',str0,repmat(str3,size(th0))),...
	       'True theta','s2','nu','rho'))
  disp(sprintf(sprintf('%s : %s ',str0,str1s),' ',th0))
end

% Don't use STRUC2ARRAY since we want them in our own order
% But see the reordering solution in OSWZEROB
if isstruct(thhats) && nargin==2
  params=thhats;
  if length(th0)>3
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str3,1,9)),...
		 'Parameters','D1','D2','g','z2','dy','dx','Ny','Nx','blurs'))
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,1,9)),...
                 ' ',params.DEL,params.g,params.z2,params.dydx,params.NyNx,params.blurs))
  else
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str3,1,5)),...
		 'Parameters','dy','dx','Ny','Nx','blurs'))
    disp(sprintf(sprintf('%s : %s ',str0,repmat(str1,1,5)),...
                 ' ',params.dydx,params.NyNx,params.blurs))
  end
else
  % Estimated values
  disp(sprintf(sprintf('%s : %s \n',str0,str1s),...
	       'Average estimated theta',mean(thhats,1)))
  % Average Hessian and Fisher matrix
   disp(sprintf('Over %i simulations, the average Hessian and the unblurred Fisher matrix are',nl))
   disp(sprintf(...
       '|%4.2f|%s apart on average (the relevant diagn file had the full information)',...
       1/100*round(100*mean(abs([avhs-Fisher]'./Fisher'*100))),'%'))
  
  % Covariance, relative, empirical, and theoretical
  disp(sprintf(sprintf('%s : %s',str0,str1s),...
	       'Observed standard deviation',std(thhats)))
  disp(sprintf(sprintf('%s : %s',str0,str1s),...
	       'Theortcl standard deviation',sqrt(diag(truecov))))
  disp(sprintf(sprintf('%s : %s',str0,str2s),...
	       'Perct of obs to pred stnddv',...
	       round(100*std(thhats)./sqrt(diag(truecov)'))))
  disp(sprintf(sprintf('%s : %s\n',str0,str2s),...
	       'Observed percent stand devn',...
	       round(100*std(thhats)./th0)))
end

% Optional output
varns={str0,str1,str2,str3};
varargout=varns(1:nargout);
