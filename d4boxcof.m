function cofs=d4boxcof
% cofs=D4BOXCOF
%
% Puts all the required coefficients for the 4-tap Daubechies wavelet
% transform on the interval in a structure array.
%
% SEE ALSO: D4BOXSTEP, D4BOXSTEPI, PRECOND4
%
% Last modified by fjsimons-at-alum.mit.edu, 03/15/2010

% Get the product filter coefficients
[h0,f0]=wc('Daubechies',2);
% Get the lowpass (H0) and highpass (H1) coefficients of the ANALYSIS
% Get the lowpass (F0) and highpass (F1) coefficients of the SYNTHESIS
[H0,H1,F0,F1]=prodco(f0,h0);

% Now the special ones for the edges, e.g. from
% http://www.pacm.princeton.edu/~ingrid/publications/54.txt
% http://www.nr.com/contrib/beylkin/D04qmf.txt
% Left, lowpass
% LLO=[          d4h1l1             d4h1l2             d4h1l3                 0                  0; ...
%                d4h2l1             d4h2l2             d4h2l3            d4h2l4             d4h2l5];
LLO=[ 0.603332511928053  0.690895531839104 -0.398312997698228                 0                  0; ...
      0.0375174604524466 0.457327659851769  0.850088102549165 0.223820356983114 -0.129222743354319];

% Left, highpass
% MINUS EVERYTHING TO COMPARE TO BEYLKIN BUT CORRECT ACCORDING TO INGRID
% LHI=[          d4g1l1             d4g1l2             d4g1l3                  0                  0; ...
%     [          d4g2l1             d4g2l2             d4g2l3             d4g2l4             d4g2l5];
LHI=[-0.796543516912183  0.546392713959015 -0.258792248333818                  0                  0; ...
      0.0100372245644139 0.122351043116799  0.227428111655837 -0.836602921223654  0.483012921773304];

% Right, lowpass
% RLO=[          d4h2r5            d4h2r4             d4h2r3             d4h2r2             d4h2r1; ...
%                     0                 0             d4h1r3             d4h1r2             d4h1r1];
RLO=[ 0.443149049637559 0.767556669298114  0.374955331645687  0.190151418429955 -0.194233407427412; ...
      0                 0                  0.230389043796969  0.434896997965703  0.870508753349866];

% Right, hipass
% MINUS EVERYTHING TO COMPARE TO BEYLKIN AND ALSO COMPARED TO INGRID
% RHI=[          d4g2r5             d4g2r4             d4g2r3             d4g2r2             d4g2r1; ...
%                     0                  0             d4g1r3             d4g1r2             d4g1r1];
RHI=[ 0.231557595006790  0.401069519430217 -0.717579999353722 -0.363906959570891  0.371718966535296;...
      0                  0                 -0.539822500731772  0.801422961990337 -0.257512919478482];

% Rearrange the coefficients for the edge calculations
% This is what will happen in D4BOXSTEP to give f(1,2,3,k+1,k+2,k+3)
LFT=[LLO ; 0 0 0 H0(4:-1:3) ; LHI ; 0 0 0 H1(4:-1:3)];

% Check the orthogonality 
difer(LFT'*LFT-eye(length(LLO)),[],[],NaN)

% This is what will happen in D4BOXSTEP to give f(k-2,k-1,k,2k-2,2k-1,2k)
RGT=[H0(2:-1:1) 0 0 0 ; RLO ; H1(2:-1:1) 0 0 0 ; RHI];

% Check the orthogonality 
difer(RGT'*RGT-eye(length(RLO)),[],[],NaN)

% Preconditioning coefficients
% On the left, forward
% HERE WE READ DOWN INGRID'S COLUMNS AND THE BEYLKIN LABELS MAKE SENSE
% LF=[         d4pc1l1                0; ...
%              d4pc1l2          d4pc2l2];
LF=[0.324894048898962  0               ; ...
    0.0371580151158803 1.00144540498130];

% On the left, inverse
% LI=[        d4ipc1l1                  0; ...
%             d4ipc1l2           d4ipc2l2];
LI=[  3.07792649138669  0                ; ...
    -0.114204567242137  0.998556681198888];

% Check the inverse 
difer(inv(LF)-LI,[],[],NaN)

% On the right, forward
% HERE THE LABELING SITS WELL WITH BEYLKIN AND WE READ INGRIDS COLUMNS WITH PERMUTATION THROUGH CENTER
% RF=[       d4pc1r1            d4pc1r2; ...
%                  0            d4pc2r2];
RF=[2.09629288435324 -0.800813234246437; ...
    0                 1.08984305289504];

% On the right, inverse
% RI=[       d4ipc1r1          d4ipc1r2; ...
%                   0          d4ipc2r2];
RI=[0.477032578540915 0.350522032550918; ...
    0                 0.917563310922261];  

% Check the inverse 
difer(inv(RF)-RI,[],[],NaN)

% Collect all of these filter coefficients into a structure
cofs.H0=H0; cofs.H1=H1; cofs.F0=F0; cofs.F1=F1; 

% Edge coefficients
cofs.LLO=LLO; cofs.LHI=LHI;
cofs.RLO=RLO; cofs.RHI=RHI; 

% Combination coefficients
cofs.LFT=LFT; 
cofs.RGT=RGT;

% Preconditioners
cofs.LF=LF;   cofs.LI=LI; 
cofs.RF=RF;   cofs.RI=RI;

