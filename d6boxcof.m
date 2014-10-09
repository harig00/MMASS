function cofs=d6boxcof
% cofs=D6BOXCOF
%
% Puts all the required coefficients for the 6-tap Daubechies wavelet
% transform on the interval in a structure array.
%
% SEE ALSO: D6BOXSTEP, D6BOXSTEPI, PRECOND6
%
% Last modified by fjsimons-at-alum.mit.edu, 03/15/2010

% disp('Obtaining coefficients')

% Get the product filter coefficients
[h0,f0]=wc('Daubechies',3);
% Get the lowpass (H0) and highpass (H1) coefficients of the ANALYSIS
% Get the lowpass (F0) and highpass (F1) coefficients of the SYNTHESIS
[H0,H1,F0,F1]=prodco(f0,h0);

% Now the special ones for the edges, e.g. from
% http://www.pacm.princeton.edu/~ingrid/publications/54.txt
% http://www.nr.com/contrib/beylkin/D06qmf.txt
% Left, lowpass
% LLO=[             d6h1l1                d6h1l2                d6h1l3               d6h1l4                     0                    0                      0                    0; ...
%                   d6h2l1                d6h2l2                d6h2l3               d6h2l4                d6h2l5               d6h2l6                      0                    0; ...
%                   d6h3l1                d6h3l2                d6h3l3               d6h3l4                d6h3l5               d6h3l6                 d6h3l7               d6h3l8]
LLO=[ 3.88899763761418e-01 -8.82078282004781e-02 -8.47841308471311e-01 3.49487436741469e-01  0                     0                     0                    0                   ; ...
     -6.21148317851952e-01  5.22527393204687e-01 -2.00008003001907e-01 3.37867348635904e-01 -3.99770770411862e-01  1.64820129734343e-01  0                    0                   ; ...
     -9.58786377870484e-03  3.71225516842464e-04  3.26009710110353e-01 8.01648164413686e-01  4.72055262027448e-01 -1.40042080963951e-01 -8.54251000994755e-02 3.52196236519723e-02]; 

% Left, highpass
% LHI=[             d6g1l1               d6g1l2                d6g1l3                d6g1l4]
%                   d6g2l1               d6g2l2                d6g2l3                d6g2l4                d6g2l5                d6g2l6                    0                     0; ...
%                   d6g3l1               d6g3l2                d6g3l3                d6g3l4                d6g3l5                d6g3l6               d6g3l7                d6g3l8; ...
LHI=[ 5.83780969578456e-01 7.93618840200530e-01  1.60955177164887e-01 -5.88417112312549e-02  0                     0                    0                     0                   ; ...
     -3.49340182426010e-01 2.98920557206519e-01 -3.28301343932461e-01 -3.32263713173123e-01  6.98249720023120e-01 -2.87878999564209e-01 0                     0                   ; ...
     -1.01505899532729e-03 3.93013301879430e-05  3.45143711309844e-02  8.48698103307436e-02 -1.33730702398544e-01 -4.60406454041029e-01 8.06893221779630e-01 -3.32671258977905e-01];

% Right, lowpass
% RLO=[            d6h3r8               d6h3r7               d6h3r6               d6h3r5               d6h3r4                d6h3r3               d6h3r2                d6h3r1; ...
%                       0                    0               d6h2r6               d6h2r5               d6h2r4                d6h2r3               d6h2r2                d6h2r1; ...
%                       0                    0                    0                    0               d6h1r4                d6h1r3               d6h1r2                d6h1r1];
RLO=[3.15093823005453e-01 7.64259199273409e-01 5.20060177812491e-01 7.70293360943509e-02 6.04255818640831e-04 -9.12473562312011e-02 -1.58758215582660e-01  8.18354184018805e-02; ...
     0                    0                    1.91450544225960e-01 4.64362767365511e-01 4.90757830674591e-01  4.96964372120369e-01  4.18999228996535e-01 -2.90407851090693e-01; ...
     0                    0                    0                    0                    5.89610106858016e-02  1.50987215328402e-01  3.82360655862306e-01  9.09684994311124e-01];

% Right, hipass
% RHI=[             d6g3r8                d6g3r7                d6g3r6                d6g3r5               d6g3r4                d6g3r3                d6g3r2                 d6g3r1; ...
%                        0                     0                d6g2r6                d6g2r5               d6g2r4                d6g2r3                d6g2r2                 d6g2r1; ...
%                        0                     0                     0                     0               d6g1r4                d6g1r3                d6g1r2                 d6g1r1]
RHI=[-1.12367571585129e-01 -2.72547235184810e-01  1.39150326725670e-01  7.59876151023844e-01  1.69441479675110e-03 -2.55869891183398e-01 -4.45179444352116e-01  2.29477548350898e-01; ...
      0                     0                     1.23073831745202e-01  2.98515239695673e-01 -7.67879574290977e-01 -9.81980008912701e-02  5.22394237872772e-01 -1.53505230705603e-01; ...
      0                     0                     0                     0                     4.07477697635819e-01 -8.04233140972164e-01 4.26562218788225e-01  -7.22194876326683e-02];

% Rearrange the coefficients for the edge calculations
% This is what will happen in D6BOXSTEP to give f(1,2,3,4,5,k+1,k+2,k+3,k+4,k+5)
LFT=[LLO ; zeros(1,4) H0(6:-1:3) ; zeros(1,6) H0(6:-1:5); LHI ; zeros(1,4) H1(6:-1:3) ; zeros(1,6) H1(6:-1:5)]; 

% Check the orthogonality 
difer(LFT'*LFT-eye(length(LLO)),[],[],NaN)

% This is what will happen in D6BOXSTEP to give f(k-4,k-3,k-2,k-1,k,2k-4,2k-3,2k-2,2k-1,2k)
RGT=[H0(2:-1:1) zeros(1,6) ; H0(4:-1:1) zeros(1,4) ; RLO ; H1(2:-1:1) zeros(1,6) ; H1(4:-1:1) zeros(1,4); RHI];

% Check the orthogonality 
difer(RGT'*RGT-eye(length(RLO)),[],[],NaN)

% Preconditioning coefficients
% On the left, forward
% HERE WE READ THE COLUMNS DOWN INGRID'S FILE AND THE LABELING SITS WELL WITH BEYLKIN
% LF=[            d6pc1l1                      0              0; ...
%                 d6pc1l2                d6pc2l2              0; ...
%                 d6pc1l3                d6pc2l3        d6pc3l3];
LF=[ 1.00794157907793e-01 0                    0               ; ...
    -5.92931024415852e-01 2.13725665320657e-01 0               ; ...            
    -1.50964976106448e-02 3.06854237778467e-02 1.00018933290722];

% On the left, inverse
% LF=[           d6ipc1l1                     0                    0; ...
%                d6ipc1l2              d6ipc2l2                    0; ...
%                d6ipc1l3              d6ipc2l3             d6ipc3l3];
LI=[ 9.92120992681742e+00  0                    0                   ; ...
     2.75240372115665e+01  4.67889524872776e+00 0                   ; ...
    -6.94679698232383e-01 -1.43546705404309e-01 9.99810702932945e-01];

% Check the inverse 
difer(inv(LF)-LI,[],[],NaN)

% On the right, forward
% HERE THE LABELING SITS WELL WITH BEYLKIN AND WE READ INGRIDS COLUMNS WITH PERMUTATION THROUGH CENTER
% RF=[           d6pc1r1               d6pc1r2               d6pc1r3; ...
%                      0               d6pc2r2               d6pc2r3; ...
%                      0                     0               d6pc3r3];
RF=[5.64221252452501e+00 -6.66336794227380e+00  2.41778315064289e+00; ...
    0                     1.73763179569470e+00 -4.65879830148413e-01; ... 
    0                     0                     1.05578252781022e+00];

% On the right, inverse
% RF=[           d6ipc1r1             d6ipc1r2              d6ipc1r3; ... 
%                       0             d6ipc2r2              d6ipc2r3; ...
%                       0                    0              d6ipc3r3];
RI=[ 1.77235436569129e-01 6.79652000611259e-01 -1.05969449845818e-01; ...
     0                    5.75495914886965e-01  2.53946179271153e-01; ...
     0                    0                     9.47164755675664e-01];

% Check the inverse 
difer(inv(RF)-RI,[],[],NaN)

% Collect all of these filter coefficients into a structure
cofs.H0=H0;   cofs.H1=H1; 
cofs.F0=F0;   cofs.F1=F1; 

% Edge coefficients
cofs.LLO=LLO; cofs.LHI=LHI;
cofs.RLO=RLO; cofs.RHI=RHI; 

% Combination coefficients
cofs.LFT=LFT; 
cofs.RGT=RGT;

% Preconditioners
cofs.LF=LF;   cofs.LI=LI; 
cofs.RF=RF;   cofs.RI=RI;
