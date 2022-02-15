function [istart,istop,dist] = finddistminima(data,signal,metric)
%FINDDISTMINIMA Finds local minima of distance metric
%   [ISTART,ISTOP,DIST] = FINDDISTMINIMA(DATA,SIGNAL,METRIC) finds local
%   minima of the specified sliding distance metric between data and
%   signal.
%   
% Copyright 2016 The MathWorks, Inc.

if size(data,2) <= size(signal,2) || ~any(strcmp(metric,{'squared','symmkl'}))
  % use direct implementation and find local minima
  [istart,istop,dist] = finddistmex(data,signal,metric);
else  
  % use correlation-based implementations
  if strcmp(metric,'symmkl')
    sdist = distsymmkl(data,signal);
  else % 'squared' (default)
    sdist = distsquared(data,signal);
  end
  
  % find all local minima (including endpoints)
  if ~isempty(sdist)
    [dist,iloc] = findpeaks(-[Inf sdist Inf]);
    istart = iloc - 1;
    istop = istart + size(signal,2) - 1;
    dist = -dist;
  else
    dist = sdist;
    istart = [];
    istop = [];
  end
end

function sdist = distsquared(data,signal)
ndata = size(data,2);
nsignal = size(signal,2);

% compute total data power E(x*x')
datapwr = movsum(sum(data.*data,1),nsignal,2,'Endpoints','discard');

% compute total signal power E(y*y')
signalpwr = sum(signal(:).*signal(:));

% compute cross power via corr: E(x*y')
crosspwr = zeros(1,ndata-nsignal+1);
for i=1:size(data,1)
  crosspwr = crosspwr + conv(data(i,:),signal(i,end:-1:1),'valid');
end

% E((x-y)*(x-y)') == E(x*x' + y*y' - x*y' - x'*y)
%  complex case can be supported via x*x' + y*y' - 2*real(x*y').
%  if rotational invariance is desired use abs() instead of real().
sdist = abs(signalpwr + datapwr - 2*crosspwr);


function sdist = distsymmkl(data,signal)
% E(x*log(x/y) + y*log(y/x)) == E((x-y)*(log(x)-log(y)) 
% Some definitions use arithmetic or geometric mean or another base for
% the logarithmic terms

ndata = size(data,2);
nsignal = size(signal,2);

logdata = log(max(data,realmin));
logsignal = log(max(signal,realmin));

% compute data terms: E(x*log(x))
dataterms = movsum(sum(data.*logdata,1),nsignal,2,'Endpoints','discard');

% compute signal terms: E(y*log(y))
signalterms = sum(signal(:).*logsignal(:));

% compute cross terms: E(x*log(y) + log(x)*y)
crossterms = zeros(1,ndata-nsignal+1);
for i=1:size(data,1)
  crossterms = crossterms + conv(data(i,:),logsignal(i,end:-1:1),'valid') ...
                          + conv(logdata(i,:),signal(i,end:-1:1),'valid');
end

% E((x-y)*(log(x)-log(y)) == E(x*log(x) + y*log(y) - x*log(y) - y*log(x)) 
sdist = abs(signalterms + dataterms - crossterms);
