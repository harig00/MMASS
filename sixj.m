function [s,C,S]=sixj(l1,l2,l3,l4,l5,l6,L,meth,C,S)
% [s,C,S]=SIXJ(l1,l2,l3,l4,l5,l6,L,meth,C,S)
%
% Wigner 6j-symbol from a database precomputed by WIGNERCYCLE.
%
% INPUT:
%
% l1,l2,l3     Top row of the Wigner 6j symbol [may be same-length vectors]
% l4,l5,l6     Bottom row of the Wigner 3j symbol [may be same-length vectors]
% L            The bandwidth of the database [default: best available]
% meth         0 Uses sparse matrices [elegant, but slow]
%              1 Performs linear search on unsorted array [slow]
%              2 Performs binary search on presorted array [default]
% C,S          The column and element vectors resulting from a previous load
%
% OUTPUT:
%
% s            The (vector of) Wigner 6j symbols
% C,S          The column and element vectors good for the next load
%
% EXAMPLE:
%
% sixj('demo1') % Should return nothing if it all works
%
% SEE ALSO: WIGNER6J, WIGNER3JM, GAUNT, WIGNER0J, GUSEINOV, ZEROJ
%
% Last modified by fjsimons-at-alum.mit.edu, 02/03/2007

% MUST want to build in conditions before even accessing database
% If NONE of the inputs satisfy the rules, that is

defval('xver',0)

if ~isstr(l1) % Not a demo
  defval('C',[])
  defval('S',[])
  
  % Method
  defval('meth',2)
  % disp(sprintf('Using method %i',meth))

  % All saved values must be integers
  if sum(mod(l1,1)) || sum(mod(l2,1)) || sum(mod(l3,1)) || ...
	sum(mod(l4,1)) || sum(mod(l5,1)) || sum(mod(l6,1))
    error('All degrees and orders must be integers for the database')
  end
  
  if isempty(C) && isempty(S)
    % Got something to do
    
    % What is the lowest of the available database bandwidths?
    Els=ls2cell(fullfile(getenv('IFILES'),'WIGNER','WIGNER6JCS-*-C'));
    for index=1:length(Els)
      EL(index)=str2num(rindeks(parse(Els{index},'-'),2));
    end
    EL=sort(EL);
    % Bandwidth of the database; keep this as low as possible
    % Need to provide for empties if not found
    fmax=find(max([l1(:)' l2(:)' l3(:)' l4(:)' l5(:)' l6(:)'])<=EL);
    if ~isempty(fmax)
      defval('L',EL(indeks(fmax,1)))
    else
      defval('L',-1)
    end
    
    % Check, identify and load database
    if any([l1(:)' l2(:)' l3(:)' l4(:)' l5(:)' l6(:)']>L)
      error('Insufficient bandwidth for database')
    end
    
    % Check rules here
    fnpl1=sprintf('%s/WIGNER6JCS-%i-C',...
		  fullfile(getenv('IFILES'),'WIGNER'),L);
    fnpl2=sprintf('%s/WIGNER6JCS-%i-S',...
		  fullfile(getenv('IFILES'),'WIGNER'),L);
    if exist(fnpl1,'file')==2 && exist(fnpl2,'file')==2
      fmt1='uint64'; fmt2='float64';
      disp(sprintf('Loading %s',fnpl1))
      C=eval(sprintf('loadb(''%s'',''%s'')',fnpl1,fmt1));
      %[C,j]=sort(C);
      % Do clear and pack for the first sort - FOR large L ONLY
      %writeb(C,fnpl1,fmt1);
      %save j
      %clear C
      %disp(sprintf('Loading %s',fnpl2))
      %S=eval(sprintf('loadb(''%s'',''%s'')',fnpl2,fmt2));
      %S=S(j);
      %writeb(S,fnpl2,fmt2); clear j
      %disp('Temporary')

      disp(sprintf('Loading %s',fnpl2))
      S=eval(sprintf('loadb(''%s'',''%s'')',fnpl2,fmt2));
      
      % Whatever happens, this better be sorted; check some entries   
      randin=unique(ceil(rand(min(100,length(C)),1)*length(C)));
      if ~all(unique(C(randin))==C(randin))
	disp('Column arrays were not properly sorted')
	[C,j]=sort(C);
	% Do clear and pack for the first sort
	writeb(C,fnpl1,fmt1);
	S=S(j);
	writeb(S,fnpl2,fmt2); clear j
	disp('Column arrays now sorted once and for all')
      end
    else
      disp('Executing WIGNERCYCLE')
      % Precompute the database
      wignercycle(L,6);
      % And have a go again
      s=sixj(l1,l2,l3,l4,l5,l6,L,meth);
    end
  else
    if isempty(L)
      error('If supplying vectors with coefficients must also specify bandwidth')
    end
    % Else have C and S from a previous load and do nothing, but check
    % There is NO check on whether the L supplied is appropriate for the
    % C & S combination that is indeed supplied
    if any([l1(:)' l2(:)' l3(:)']>L)
      error('Insufficient bandwidth for database')
    end
  end
  
  % Find running index into this matrix
  % There must be the same number of all of the input degrees
  for ix=1:length(l1)
    rind(ix)=1+l1(ix)*(L+1)^0+...
	     l2(ix)*(L+1)^1+...
	     l3(ix)*(L+1)^2+...
	     l4(ix)*(L+1)^3+...
	     l5(ix)*(L+1)^4+...
	     l6(ix)*(L+1)^5;
  end
  % Initialize results vector
  s=repmat(NaN,1,length(rind));

  switch meth
   case 0
    % Turn it into the properly indexed sparse matrix
    % This step takes some of time but most of all, memory
    W=sparse(1,C,S,1,(L+1)^6);

    % Extract the Wigner 6j-symbol
    s=W(1,rind);
   case 1
    for index=1:length(rind)
      posi=find(C==rind(index));
      if ~isempty(posi)
	s(index)=S(posi);
      else 
	s(index)=0;
      end
    end
   case 2
    % Binary search algorithm on sorted arrays
    for index=1:length(rind)
      posi=binsearch(C,rind(index));
      if ~isempty(posi)
	s(index)=S(posi);
      else 
	s(index)=0;
      end
    end
   otherwise
    error('Specify valid method')
  end
  % Check some special cases by Edmonds' rules
  if xver==1 && length(l1)==1 && length(l2)==1 && length(l3)==1
    es=l1+l2+l3;
    ex=l2*(l2+1)+l3*(l3+1)-l1*(l1+1);
    if l4==0 & s~=0
      difer(s-(-1)^es/sqrt((2*l2+1)*(2*l3+1)))
    elseif l4==1 & s~=0
      difer(s-2*(-1)^(es+1)*ex...
	    /sqrt(2*l2*(2*l2+1)*(2*l2+2)*2*l3*(2*l3+1)*(2*l3+2)))
    elseif l4==2 & s~=0
      difer(s-(2*(-1)^es*(3*ex*(ex-1)-4*l2*(l2+1)*l3*(l3+1)))...
	    /sqrt((2*l2-1)*2*l2*(2*l2+1)*(2*l2+2)*(2*l2+3)*...
		  (2*l3-1)*2*l3*(2*l3+1)*(2*l3+2)*(2*l3+3)))
    end
    disp('Passed excessive verification of Edmonds'' special cases')
  end
elseif strcmp(l1,'demo1')
  l1=0:15;
  l2=8; l3=7; l4=13; l5=15; l6=15;
  W=wigner6j(l2,l3,l4,l5,l6);
  l2=repmat(l2,size(l1));
  l3=repmat(l3,size(l1));
  l4=repmat(l4,size(l1));
  l5=repmat(l5,size(l1));
  l6=repmat(l6,size(l1));
  difer(W-sixj(l1,l2,l3,l4,l5,l6))
  difer(W-[0 (1/8)*sqrt(259/5270) -(11/24)*sqrt(407/152830) ...
	   (1/204)*sqrt(54131/116870) (1/68)*sqrt(23199/23374) ... 
	   -(1/1768)*sqrt(142709/186) (59/8840)*sqrt(209/5394) ...
	   (997/2210)*sqrt(33/41354) -(1429/3910)*sqrt(77/70122) ...
	   -(7747/15640)*sqrt(1/70122) (269/15640)*sqrt(3157/5394)...
	   -(83/70380)*sqrt(13981/290) -(1/70380)*sqrt(88398583/310) ...
	   (149/1360680)*sqrt(4790071/310) (13/240120)*sqrt(368467/1054) ...
	   -(13/40455)*sqrt(2579269/782)]) 
  % Try various methods as well
  l1=0:6;
  l2=8; l3=7; l4=6; l5=4; l6=2;
  W=wigner6j(l2,l3,l4,l5,l6);
  l2=repmat(l2,size(l1));
  l3=repmat(l3,size(l1));
  l4=repmat(l4,size(l1));
  l5=repmat(l5,size(l1));
  l6=repmat(l6,size(l1));
  difer(W-[0 0 -(1/42)*sqrt(11/34) (1/21)*sqrt(209/442) -(5/21)*sqrt(19/442) ...
	   (5/6)*sqrt(19/4862) -(1/78)*sqrt(133/17)])
  difer(W-sixj(l1,l2,l3,l4,l5,l6))
  difer(sixj(4,0,4,13,9,13)-1/9/sqrt(3))
  difer(sixj(30,30,30,30,30,30)-4.102353215741e-04)
end

