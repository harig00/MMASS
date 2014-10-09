function varargout=covthos(th,params,k,scl)
% [covth,covlin,F,Flin]=covthos(th,params,k,scl)
%
% Calculates the entries of the covariance matrix of the estimate, indeed
% the inverse Fisher matrix, hopefully close to the expectation of the
% Hessian of the actual simulations. As seen in Olhede & Simons for the
% UNCORRELATED Forsyth loading model.
%
% INPUT:
%
% th       The five-parameter vector (true or estimated) [scaled]:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% k        Wavenumber(s) at which the Fisher matrix is evaluated [1/m]
% scl      The vector with any scalings applied to the parameter vector
%
% OUTPUT:
%
% covth    The covariance matrix between the parameters
% covlin   The unwrapped covariance matrix, as printed
% F        The SCALED Fisher matrix
% Flin     The unwrapped SCALED Fisher matrix, as printed
%
% Last modified by fjsimons-at-alum.mit.edu, 09/26/2012

% Default scaling is none
defval('scl',ones(size(th)))

% Scale up the parameter vector for the proper likelihood and score
th=th.*scl;

% First, the Fisher matrix at each wavenumber, unwrapped
mcF=Fisherkos(k,th,params);

% Take the expectation and put the elements in the right place
for ind=1:length(mcF)
  mcF{ind}=nanmean(mcF{ind});
end

% The full Fisher matrix
% These will become the variances
F(1,1)=mcF{1};
F(2,2)=mcF{2};
F(3,3)=mcF{3};
F(4,4)=mcF{4};
F(5,5)=mcF{5};

% These will be the covariances of D with others
F(1,2)=mcF{6};
F(1,3)=mcF{7};
F(1,4)=mcF{8};
F(1,5)=mcF{9};

% These will be the covariances of f2 with others
F(2,3)=mcF{10};
F(2,4)=mcF{11};
F(2,5)=mcF{12};

% These will be the covariances of s2 with others
F(3,4)=mcF{13};  
F(3,5)=mcF{14};

% These will be the covariances of nu with others
F(4,5)=mcF{15};

% Construct the UNSCALED Fisher matrix
F=[triu(F)'+triu(F)-diag(diag(F))]; 

% Any NaN's above came from zeroes in the wavenumbers
warning off MATLAB:nearlySingularMatrix
% Note that only the half plane is statistically independent so when
% working in the whole plane our variance would have looked too small
isfullplane=1; 
ishalfplane=0;
if isfullplane==1
  disp('I am assuming that your wavenumbers are the entire plane')
  covth=inv(F)/length(k(~~k))*2;
elseif ishalfplane==1
  disp('I am assuming that your wavenumbers are the half plane')
  covth=inv(F)/length(k(~~k));
end
warning on MATLAB:nearlySingularMatrix

% This as printed for reading by DIAGNOS
% Unwrapped and in the unscaled dimensions of the problem
covlin=tril(covth); covlin=covlin(~~covlin)';
% Unwrapped and in the scaled dimensions of the problem
F=[scl(:)*scl(:)'].*F;
Flin=tril(F); Flin=Flin(~~Flin)';

% Output
varns={covth,covlin,F,Flin};
varargout=varns(1:nargout);
