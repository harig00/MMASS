function [ps,yl]=succaplot(tim,xs,diro,scl,dirmi)
% [ps,yl]=SUCCAPLOT(tim,xs,diro,scl,dirmi)
%
% Plots the cumulative sum of a matrix, 
% appropriately offset for upgoing or
% downgoing plots.
%
% 'diro' 'up' or 'dn'
% 'scl'  scale factor for offset [default: 1]
% 'dirmi' minimum offset [default: half the largest]
%
% ONLY TESTED FOR 'UP' THUS FAR
%
% Last modified by fjsimons-at-alum.mit.edu, March 25th, 2003

defval('diro','up')
defval('scl',1)

% Make sum and demean
cxs=cumsum(xs,2);
cxs=cxs-mean(cxs(:));
mins=min(cxs);
maxs=max(cxs);

% Calculate offsets
offs(1)=0;
for index=1:size(xs,2)
  if index<size(xs,2)
    if diro=='up'
      offs(index+1)=offs(index)+maxs(index)-mins(index+1);
    elseif diro=='dn'
      offs(index+1)=ofss(index)+mins(index)-maxs(index+1);
    end
  end
end
offs=offs*scl;
doffs=diff(offs);

defval('dirmi',max(doffs)/2)

% Plots signal
for index=1:size(xs,2)
  if index>1
    if doffs(index-1)<dirmi
      offs(index:end)=offs(index:end)+(dirmi-doffs(index-1));
    end
  end
  yl(index)=offs(index);
  ps(index)=plot(tim,cxs(:,index)+offs(index));
  hold on
end
hold off
set(gca,'Ytick',offs,'YTickL',fliplr([1:size(xs,2)]))
ylim([offs(1)-dirmi offs(end)+dirmi])
