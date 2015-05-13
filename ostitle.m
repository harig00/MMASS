function t=ostitle(ah,params,ovall,numsim)
% t=OSTITLE(ah,params,ovall,numsim)
% 
% Makes a decent overview graph title
%
% INPUT:
%
% ah          Axis handles over which to put the title
% params      Parameter structure of the experiment
% ovall       Overall plot title for identification
% numsim      Number of simulations to work in the title
%
% OUTPUT:
%
% t           Title handle
%
% Last modified by fjsimons-at-alum.mit.edu, 10/24/2014

defval('ovall','')
defval('numsim','')

if numsim>1
  mlest=sprintf('%i MLE simulations ',numsim);
else
  mlest=sprintf('%i MLE simulation ',numsim);
end

if ~isempty(ovall)
  str9='%s\n\n%s %s_1 = %i ; %s_2 = %i %s ; z_2 = %i km ; %ix%i grid ; %ix%i km ; blur %i';
  str5='%s\n\n%s ; %ix%i grid ; %ix%i km ; blur %i';
else
  str9='%s%swith %s_1 = %i ; %s_2 = %i %s\n z_2 = %i km ; %ix%i grid ; %ix%i km ; blur %i';
  str5='%s%swith %ix%i grid ; %ix%i km ; blur %i';
end

% Move away from numbered to named access to the params...
struct2var(params)

% Conditions pertain to cases with and without the kiso parameter
% present, which really should be recorded in the simulations
if length(fieldnames(params))>=7
  t=supertit(ah,sprintf(str9,ovall,mlest,...
      '\Delta',DEL(1),'\Delta',DEL(2),'kg m^{-3}',...
      z2/1e3,NyNx,dydx.*NyNx/1e3,blurs));
elseif length(fieldnames(params))>=4
  t=supertit(ah,sprintf(str5,ovall,mlest,...
      NyNx,dydx.*NyNx/1e3,blurs));
else
  error('OSTITLE: Not the expected number of fixed parameters!')
end




