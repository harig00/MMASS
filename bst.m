function [qual,sncal]=bst(D,T,T0,T1)
% [qual,sncal]=bst(D,T,T0,T1)
%
% Comes up with a quality criterion based on wavelet coefficients.
% Now check out the means of the scalograms between the triggering points
% Signal-to-noise ratio is defined by the ratio of the power
% between T0 and T1 at scales 2, 3, 4 compared to the time interval
% before it at the same scales.
%
% INPUT:
%
% D       Cell array with wavelet coefficients
% T       Cell array with time axis
% T0      Triggered time (s)
% T1      Detriggered time (s)
%
% OUTPUT:
%
% qual   0 Definitely not a good selection
%        1 Clear waveform, low noise
%        2 Might be a good selection, but not great
% sncal  Estimate of the signal-to-noise ratio
%
% SEE ALSO: BFT.
%
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2010

for j=1:length(D)
  b(j)=mean(abs(D{j}(T{j}>T0 & T{j}<T1)));
end

% Lots of power in scales 3 and 4
if all(b(3)>b([1 2 5])) & all(b(4)>b([1 2 5]))
  qual=1;
  % Much more in scales 1, 2 and 3 than in 4 or 5
elseif all(b(1)>b([4 5])) & all(b(2)>b([4 5])) &  all(b(3)>b([4 5]))
  qual=0;
else
  qual=2;
  % Maybe some power at scale 5 too; but otherwise 3 and 4 win out
  if all(b(3)>b([1 2])) & all(b(4)>b([1 2]))
    qual=1;
  end
  % Lots of power in scales 2 and 3
  if all(b(2)>b([1 4 5])) & all(b(3)>b([1 4 5]))
    qual=1;
  end
end
% If there is less power in the promising interval than in the interval
% before it drop the thing altogether
for j=1:length(D)
  b2(j)=mean(abs(D{j}(T{j}>T0-(T1-T0) & T{j}<T0)));
end
if any(b(3:4)<b2(3:4))
  qual=0;
end
% Get rid of short T-phases?
% This didn't help one bit
% if all(b>8*b2); qual=0; end

% Estimate of signal-to-noise ratio
sncal=mean(b([2 3 4]))./mean([b2([1:5])]);
