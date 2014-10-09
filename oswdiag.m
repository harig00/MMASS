function oswdiag(fid3,fmt1,fmt3,logli,gr,hs,thhat,thini,scl,ts,e,o,vHxs)
% OSWDIAG(fid3,fmt1,fmt3,logli,gr,hs,thhat,thini,scl,ts,e,o,vHxs)
%
% Writes an entry in the DIAGN file for Olhede & Simons (2013) suite
%
% Last modified by fjsimons-at-alum.mit.edu, 10/02/2014

% There may be one, there may be two samples variances incoming
fmtx=repmat(' %9.3e',1,length(vHxs));

% Turn the Hessian into something good
hessi=tril(full(hs)); hessi=hessi(~~hessi);

% Now also stick in the spatial variance estimated poorly
% This whole thing gets mangled but I like to keep 'fmt1' clean
% Borrow the format from the format of the sigma^2 parameter
% First line: The initial guess, scaled back to proper units
fprintf(fid3, fmt1,                     thini.*scl      );
% Second line: The estimate, scaled back to proper units
fprintf(fid3,[fmt1(1:end-2) fmtx '\n'],[thhat.*scl vHxs]);
% Third line: Three diagnostics, likelihood, diagnostic, scales, and derivatives
fprintf(fid3,fmt3,...
	round(ts),e,o.iterations,...
	logli,o.firstorderopt,...
	scl,gr,hessi);
% Fourth line: Empty (only visuals, no effect on reading)
fprintf(fid3,'\n');
