function [newmat,trimi]=trimit(oldmat,perc,maxi,col)
% [newmat,trimi]=TRIMIT(oldmat,perc,maxi,col)
%
% Trims rows of a data set based on a percentage criterion applied to the
% values in a certain column.
%
% INPUT:
%
% oldmat      Input data matrix, observations are rows, columns are values
% perc        Percentage for trimming decision [default: 100, keep all]
% maxi        0 keep the central perc% portion of the data
%             +1 keep the lowest perc% of the data (trim off right)
%             -1 keep the highest perc% of the data (trim off left)
% col         Column indices for trimming decision [default: all]
%             These are considered independently, recursively
%
% OUTPUT:
%
% newmat      The trimmed results, of smaller size [if vector, perc% remain]
% trimi       The indices of the surviving rows; logical or indicial
%
% SEE ALSO:
% 
% WINSOR
%
% Last modified by fjsimons-at-alum.mit.edu, 02/06/2015

defval('perc',100)
defval('maxi',0)
defval('col',1:size(oldmat,2))

if length(col)>1
  trimo=[1:size(oldmat,1)]';
  for index=1:length(col)
    [oldmat,trimi]=trimit(oldmat,perc,maxi,col(index));
    trimo=trimo(trimi);
  end
  newmat=oldmat;
  % A numeric row index
  trimi=trimo;
else
  switch maxi
   case 0
    perx=50+[-1 1]*perc/2;
   case +1
    perx=[0 perc];
   case -1
    perx=[100-perc 100];
  end
  if perc==0
    [newmat,trimi]=deal([]);
  elseif perc==100
    newmat=oldmat;
    trimi=logical(ones(1,size(oldmat,1)));
  else
    cutpts=prctile(oldmat(:,col),perx);
    % A logical row index
    trimi=[oldmat(:,col)>=cutpts(1)] & [oldmat(:,col)<=cutpts(2)];
    newmat=oldmat(trimi,:);
  end
  % Report on the trimmings
  disp(sprintf('%i%s trimmed',round(length(newmat(:))/length(oldmat(:))*100),'%'))
end
