function oswzerob(fid,th0,params,options,bounds,fmt1,fmt2)
% OSWZEROB(fid,th0,params,options,bounds,fmt1,fmt2)
% 
% Writes the beginning of a THZRO file as we have come to use it 
%
% INPUT:
%
% fid        The file id (from OSOPEN)
% th0        The parameter vector (see, e.g., SIMULOS)
% params     A structure with the known constants (see, e.g. SIMULOS)
% options    The options used by the optimization procedure
% bounds     The bounds used by the optimization procedure
% fmt1...    Strings containing formatting instructions (from OSOPEN)
%
% SEE ALSO: 
%
% OSWZEROE, OSRZERO
%
% Last modified by fjsimons-at-alum.mit.edu, 02/09/2015

% Commit the truth to file
fprintf(fid,'%s\n','the true parameter vector');
fprintf(fid,fmt1,th0);

% Commit the parameters of the experiment to file
fprintf(fid,'%s\n','the fixed experimental parameters');

%fprintf(fid,fmt2,struct2array(params));
% Rather, these need to be ordered to yield to the format
fulls={'DEL','g','z2','dydx','NyNx','blurs','kiso'};
[~,i]=ismember(fulls,fieldnames(params));
jk=struct2cell(params);
fprintf(fid,fmt2,[jk{i(~~i)}]);

% Convert the bounds to something printable
fprintf(fid,'%s\n','the bounds, if any');
if ~isempty(bounds)
  struct2str(cell2struct(bounds,...
		    {'A',  'B'  ,... % Linear Inequalities
		    'Aeq','Beq',... % Linear Equalities
		    'LB',...        % Lower Bounds
		    'UB',...        % Upper Bounds
		    'NONLCON'},...   % Nonlinear Inequalities
			 2),fid);
else
  % We need at least one colon on the next line
  struct2str(cell2struct({'None'},{'Bounds'},2),fid)
end

% Commit the parameters of the optimization to file
fprintf(fid,'%s\n','the optimization options');
struct2str(options,fid)
