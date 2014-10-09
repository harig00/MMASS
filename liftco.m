function [Pc,Uc,KPc,KUc,K,U,P]=liftco(tipe,N,M)
% [Pc,Uc,KPc,KUc,K,U,P]=LIFTCO(tipe,N,M)
%
% Like WC except also calculates Laurent polynomial matrices.
% Check with LSINFO...
%
% Comes up with Prediction and Update coefficients
% for a lifting implementation of the wavelet transform
% of a certain type and with parameters N and M
% where N is the vanishing moments of the primal wavelet
%       N-1  is the order of polynomial cancellation
%       M is the vanishing moments of the dual wavelet
%       M-1  is the order of preserved moments
% These coefficients assume uniform sample spacing and perform a two-phase
% "Lazy" preprocessing step. 
% Also returns the appropriate K,U,P which should yield
% the analysis polyphase matrix upon Laurent multiplication.
%
% CDF is the Cohen-Daubechies-Feauveau construction
% CDF(1,1) is the Orthogonal Haar bank
% CDF(4,2) is the cubic B-spline
%
% See also MAKEWC and WC.
%
% Last modified by fjsimons-at-alum.mit.edu, 18.04.2006

defval('tipe','CDF')
defval('N',1)
defval('M',1)
 
% Collect the operators required
% They might be cell arrays themselves
[h0,f0,Pc,Uc,KPc,KUc]=wc(tipe,[N M]);

% Verify factorization with a Laurent polynomial
% multiplication.
if ~iscell(Pc) & ~iscell(Uc) 
  K=[{KPc} {0} ; {0} {KUc}];
  U=[{1} {Uc} ; {0} {1}];
  P=[{0} {1} ;  {1} {-Pc}]; 
else
  [K,U,P]=deal([]);
end



