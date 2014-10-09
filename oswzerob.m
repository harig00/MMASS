function oswzerob(fid,th0,params,options,fmt1,fmt2)
% OSWZEROB(fid,th0,params,options,fmt1,fmt2)
% 
% Writes the beginning of a THZRO file as we have come to use it 
%
% INPUT:
%
% fid        The file id (from OSOPEN)
% th0        The parameter vector (see, e.g., SIMULOS)
% params     A structure with the known constants (see, e.g. SIMULOS)
% options    The options used by the optimization procedure
% fmt1...    Strings containing formatting instructions (from OSOPEN)
%
% SEE ALSO: 
%
% OSWZEROE, OSRZERO
%
% Last modified by fjsimons-at-alum.mit.edu, 10/06/2014

% Commit the truth to file
fprintf(fid,'%s\n','the true parameter vector');
fprintf(fid,fmt1,th0);

% Commit the parameters of the experiment to file
fprintf(fid,'%s\n','the fixed experimental parameters');
fprintf(fid,fmt2,struct2array(params));

% Commit the parameters of the optimization to file
fprintf(fid,'%s\n','the optimization options');
struct2str(options,fid)
