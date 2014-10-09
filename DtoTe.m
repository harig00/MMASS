function Te=DtoTe(D,E,nu)
% Te=DTOTE(D,E,nu)
%
% Converts flexural rigidity to elastic thickness in SI. 
%
% INPUT:
%
% D    Flexural rigidity [Nm]
% E    Young's modulus [Pa] (defaulted)
% nu   Poisson's ratio (defaulted)
%
% OUTPUT:
%
% Te   Effective elastic thickness [m]
%
% Last modified by fjsimons-at-alum.mit.edu, 08/12/2013

% Young's modulus [Pa]
defval('E',1.4e11);
% Poisson's ratio
defval('nu',0.25);

% Effective elastic thickness [m]
Te=[D*12*(1-nu^2)/E].^(1/3);
