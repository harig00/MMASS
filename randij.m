function siz=randij(nup)
% siz=RANDIJ(nup)
%
% INPUT:
%
% nup       A certain total number of elements
% 
% OUTPUT:
%
% siz       A nontrivial random size of a matrix with nup elements

defval('nup',180); 
siz=shuffle(factor(nup)'); 
rnp=randi(length(siz)-1); 
siz=[prod(siz(1:rnp)) prod(siz(rnp+1:end))];


