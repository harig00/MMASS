function A=realize(A,tolp)
% A=realize(A,tolp)
%
% Gets rid of small imaginary components of a matrix and reports on the
% committed error before doing it, within tolerance.
%
% INPUT:
%
% A      The matrix in question
% top    The tolerance in percentage [default: 0.001]
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2014

defval('tolp',0.001);

% Watch the case where it's just pi

% Which function is invoking us
ds=dbstack;

% Report on small imaginary parts if any
cpxity=100*nanmean(abs(imag(A(:))))./nanmean(abs(real(A(:))));
if cpxity<tolp
  A=real(A);
else
  % Still do it, but report on it!
  A=real(A);
  disp(sprintf('%s Imaginary/real parts ratio is %5.0e%s',...
	       upper(ds(2).name),cpxity,'%'))
end
