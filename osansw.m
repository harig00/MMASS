function [answ,answs]=osansw(th0,truecov,E,v)
% [answ,answs]=OSANSW(th0,truecov,E,v)
%
% Constructs some reporting strings for MLEOS, MLEROS, MLEROS, MLEOSL, 
% thus for correlated and uncorrelated loading scenarios as well as the
% single-field estimate
%
% INPUT:
%
% th0        True parameter vector
% truecov    The theoretical parameter covariance
% E          Young's modulus [Pa]
% v          Poisson's ratio
%
% OUTPUT:
%
% answ       A cell with the "answers" if you will
% answs      Appropriate formatting strings for a plot title
%
% Last modified by fjsimons-at-alum.mit.edu, 08/30/2013

defval('E',[])
defval('v',[])

% The plus-minus sign
pm=str2mat(177);

% Dual-field analysis with no business for E or v
if ~isempty(E) && ~isempty(v)
  % Conversion to Te [m]
  Te=DtoTe(th0(1),E,v);
  
  % This will come awfully close but MUST PUT IN THE MEAN!
  stdTe=std(DtoTe(th0(1)+randn(10000,1)*sqrt(truecov(1,1)),E,v));
  stdf2=sqrt(truecov(2,2));
  
  % Here is the "delta method", which is approximate
  stdTedelta=sqrt(truecov(1,1))/th0(1)^(2/3)/3*[12*(1-v^2)/E]^(1/3);

  disp(sprintf('Stdv of the Te predicted by simulation:  %3.2f km',stdTe/1000))
  disp(sprintf('Stdv of the Te predicted by delta meth:  %3.2f km\n',stdTedelta/1000))
  
  % True parameters and their theoretical uncertainties
  answ{1}=sprintf('E = %g',E);
  answ{2}=sprintf('%s = %g','\nu',v);
  answ{3}=sprintf('T_e = %4.1f %s %3.1f km',...
		  round(Te/1e3*10)/10,pm,round(stdTe/1e3*10)/10);
  answ{4}=sprintf('f^2 = %g %s %g',...
		  th0(2),pm,round(stdf2*1e3)/1e3);
end
  
% Some common strings
s2s='%s^2 = %6.4f %s %6.4f';
nus='%s = %5.3f %s %5.3f';
rhs='%s = %i %s %i';

if length(th0)==3
  % Single-field analysis
  stds2=sqrt(truecov(1,1));
  stdnu=sqrt(truecov(2,2));
  stdrho=sqrt(truecov(3,3));

  answ{1}=sprintf(s2s,'\sigma',th0(1),pm,stds2);
  answ{2}=sprintf(nus,'\nu',   th0(2),pm,stdnu);
  answ{3}=sprintf(rhs,'\rho',  round(th0(3)),pm,round(stdrho));

  answs=['%s ; %s ; %s'];
elseif length(th0)==5
  % Uncorrelated analysis
  stds2=sqrt(truecov(3,3));
  stdnu=sqrt(truecov(4,4));
  stdrho=sqrt(truecov(5,5));
  
  answ{5}=sprintf(s2s,'\sigma',th0(3),pm,stds2);
  answ{6}=sprintf(nus,'\nu',   th0(4),pm,stdnu);
  answ{7}=sprintf(rhs,'\rho',  round(th0(5)),pm,round(stdrho));

  answs=['%s ; %s ; %s ; %s\n %s ; %s ; %s'];
elseif length(th0)==6
  % Correlated analysis
  stdr=sqrt(truecov(3,3));
  stds2=sqrt(truecov(4,4));
  stdnu=sqrt(truecov(5,5));
  stdrho=sqrt(truecov(6,6));
  
  answ{5}=sprintf('r = %g %s %g',	th0(3),pm,round(stdr*1e3)/1e3);
  answ{6}=sprintf(s2s,'\sigma',	th0(4),pm,stds2);
  answ{7}=sprintf(nus,'\nu',	th0(5),pm,stdnu);
  answ{8}=sprintf(rhs,'\rho',	round(th0(6)),pm,round(stdrho));
  
  answs=['%s ; %s ; %s ; %s ; %s\n %s ; %s ; %s'];
end
