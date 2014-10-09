function t=ostitle(ah,params,ovall,numsim)
% t=OSTITLE(ah,params,ovall,numsim)
% 
% Makes a decent overview graph title
%
% INPUT:
%
% ah          Axis handles over which to put the title
% params      Parameters that make it into the title string, vector or structure
% ovall       Overall plot title for identification
% numsim      Number of simulations to work in the title
%
% OUTPUT:
%
% t           Title handle
%
% Last modified by fjsimons-at-alum.mit.edu, 07/02/2014

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

% This is a bit too dependent on the order they come out
if isstruct(params)
  params=struct2array(params);
end

if isfield(params,'kiso')
  % Reorganize the kiso parameter to be at the end if it exists
  params=[params(2:end) params(1)];
end

% Should move away from numbered to named access to the params...
% Conditions pertain to cases with and without the kiso parameter
% present, which really should be recorded in the simulations
if length(params)>=9
  t=supertit(ah,sprintf(str9,ovall,mlest,...
      '\Delta',params(1),'\Delta',params(2),'kg m^{-3}',...
      params(4)/1e3,params(7:8),round(params(7:8).*params(5:6)/1e3),params(9)));
elseif length(params)>=5
  t=supertit(ah,sprintf(str5,ovall,mlest,...
      params(3:4),round(params(3:4).*params(1:2)/1e3),params(5)));
else
  error('OSTITLE: Not the expected number of fixed parameters!')
end




