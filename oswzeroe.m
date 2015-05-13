function oswzeroe(fid,scl,avH,good,F,covF,fmtc,fmte,fmtf)
% OSWZEROE(fid,scl,avH,good,Flin,covlin,fmtc,fmte,fmtf)
% 
% Writes the (e)nd of a THZRO file as we have come to use it 
%
% INPUT:
%
% fid          The file id (from FOPEN)
% scl          The scaling factors for the Fisher matrix 
% avH          The full-form sample-average Hessian matrix
% good         The number of samples over which this was averaged
% F            The full-form scaled Fisher matrix
% covF         The full-form theoretical covariance matrix, basd on F 
% fmt{c,e,f}   Strings containing formatting instructions    
%
% SEE ALSO: 
%
% OSWZEROB, OSRZERO
%
% Last modified by fjsimons-at-alum.mit.edu, 02/09/2015

% Print the scaling of the theoretical values
fprintf(fid,'%s\n','the scaling factors');
fprintf(fid,fmtc,scl);

% Now print the theoretical covariance to file also
fprintf(fid,'%s\n','the theoretical covariance');
fprintf(fid,fmtf,trilos(covF));

% Print the scaled Fisher matrix, the expected value of the Hessian
fprintf(fid,'%s\n','the unblurred scaled Fisher matrix');
fprintf(fid,fmte,trilos(F));

% Print the observed average of the Hessians (over last set of runs)
fprintf(fid,'%s\n',...
        sprintf('the scaled Hessian matrix averaged over whichever %i runs last output',...
                good));
fprintf(fid,fmte,trilos(avH));
