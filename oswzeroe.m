function oswzeroe(fid,scl,avlin,good,Flin,covlin,fmtc,fmte,fmtf)
% OSWZEROE(fid,scl,avlin,good,Flin,covlin,fmtc,fmte,fmtf)
% 
% Writes the (e)nd of a THZRO file as we have come to use it 
%
% INPUT:
%
% fid          The file id (from FOPEN)
% scl          The scaling factors for the parameter vector
% avlin        The sample average Hessian matrix
% good         The number of samples over which this was averaged
% Flin         The expectation of the Hessian matrix (the Fisher matrix)
% covlin       The theoretical covariance matrix
% fmt{c,e,f}   Strings containing formatting instructions    
%
% SEE ALSO: 
%
% OSWZEROB, OSRZERO
%
% Last modified by fjsimons-at-alum.mit.edu, 10/06/2014

% Print the single scaling of the theoretical values that was used 
fprintf(fid,'%s\n','the experimental scaling factors');
fprintf(fid,fmtc,scl);

% Now print the theoretical covariance to file also
fprintf(fid,'%s\n','the theoretical covariance');
fprintf(fid,fmtf,covlin);

% Print the theoretical average of the Hessians
fprintf(fid,'%s\n','the scaled Fisher matrix');
fprintf(fid,fmte,Flin);

% Print the observed average of the Hessians (over last set of runs)
fprintf(fid,'%s\n',...
        sprintf('the scaled Hessian matrix averaged over whichever %i runs last output',...
                good));
fprintf(fid,fmte,avlin);
