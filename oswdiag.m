function oswdiag(fid3,fmt1,fmt3,logli,gr,hs,thhat,thini,scl,ts,e,o,vHxs,momx,covh)
% OSWDIAG(fid3,fmt1,fmt3,logli,gr,hs,thhat,thini,scl,ts,e,o,vHxs,momx,covh)
%
% INPUT:
%
% fid3,fmt1,fmt3   Format strings
% K                Number of independent wavenumbers, equal to 
%                  length(k(~~k))/2 [entire plane]
% logli,gr,hs      Likelihood, score, Hessian
% thhat,thini,scl  Estimates, initial values, scales
% ts,e,o           Optimization timing, exit flag, optimality
% vHxs             Spatial (sample) variance
% momx             Moments of the quadratic portion of the likelihood
% covh             The asymptotic covariance matrix for the estimate
%
% Writes an entry in the DIAGN file for the Olhede & Simons (2013) suite
%
% Last modified by fjsimons-at-alum.mit.edu, 10/25/2014

% There may be one, there may be two sample variances incoming
fmtx=repmat(' %9.3e',1,length(vHxs));

% First line: The initial guess, scaled back to proper units
fprintf(fid3, fmt1,                     thini.*scl      );
% Second line: The estimate, scaled back to proper units, the sample variance
fprintf(fid3,[fmt1(1:end-2) fmtx '\n'],[thhat.*scl vHxs]);
% Other lines: Three diagnostics, likelihood, diagnostic, moments, scales, and derivatives
fprintf(fid3,fmt3,...
	round(ts),e,o.iterations,...
	logli,o.firstorderopt,momx,...
	scl,gr,trilos(hs),trilos(covh));
% Last line: Empty (only visuals, no effect on reading)
fprintf(fid3,'\n');
