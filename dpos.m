function dprod=dpos(DEL,varargin)
% dprod=dpos(DEL,POW1,POW2,...)
%
% Forms the product of the powers of the elements of a vector involving
% the densities as needed in Olhede & Simons.
%
% INPUT:
%
% DEL        A 1xN vector (in casu, of two density contrasts)
% POW1       The first power
% POW2       The second power, ..., until POWN     
%
% OUTPUT:
%
% dprod      DEL(1)^POW1*DEL(2)^POW2*...*DEL(N)^POWN;
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2012

POW=cat(2,varargin{:});
dprod=prod(DEL.^POW);
