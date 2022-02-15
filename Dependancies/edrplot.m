function edrplot(x, y, ix, iy, dist)
%EDRPLOT Helper function for plotting dynamically time warped signals
%
%   For vector X and Y, EDRPLOT plots the original signals and aligned
%   signals, where the original signals are placed in a (smaller) plot.
%
%   For matrix X and Y, EDRPLOT plots X and Y as images (via IMAGESC).  
%   To facilitate visual alignment of features, EDRPLOT additionally
%   plots X and Y together with 50% transparency where they overlap.

%   Copyright 2015 The MathWorks, Inc.

if isvector(x) && isvector(y)
  plotVectorAlignment(x, y, ix, iy, dist);
else
  plotMatrixAlignment(x, y, ix, iy, dist);
end

%-------------------------------------------------------------------------
function plotVectorAlignment(x, y, ix, iy, dist)
% plot original signals in the first plot
if isreal(x)
  subplot('Position',[0.13 0.74 0.78 0.17]); 
  plot(1:numel(x),x,'.-', ...
       1:numel(y),y,'.-')
else
  subplot('Position',[0.13 .8 .75 .15]);
  plot3(1:numel(x),real(x),imag(x),'.-', ...
        1:numel(y),real(y),imag(y),'.-');
  view(-1.9,19.8);
  ylabel('real');
  zlabel('imag');
end
title(getString(message('signal:edr:OriginalSignals')));

% plot aligned signals in the second plot
if isreal(x)
  subplot('Position',[0.13 0.11 0.78 0.49]);
  plot(1:numel(ix),x(ix),'.-', ...
       1:numel(iy),y(iy),'.-');
else
  subplot('Position',[0.13 0.11 .75 .55]);
  plot3(1:numel(ix),real(x(ix)),imag(x(ix)),'.-', ...
        1:numel(iy),real(y(iy)),imag(y(iy)),'.-');
  view(-3.5,13.4);
  ylabel('real');
  zlabel('imag');
end
titleStr = sprintf('%s (%s: %i)', ...
   getString(message('signal:edr:AlignedSignals')), ...
   getString(message('signal:edr:Edit')), dist);
title(titleStr);

%-------------------------------------------------------------------------
function plotMatrixAlignment(x, y, ix, iy, dist)
if ~isreal(x) || ~isreal(y)
  % split each matrix into its real and imaginary portions
  x = [real(x); imag(x)]; y = [real(y); imag(y)];
end

% Plot original signals in first column of plots
subplot('Position',[0.07 0.71 0.4 0.19])
imagesc(x);
title(getString(message('signal:edr:OriginalSignal','X')));
xlim([1 max(size(x,2),size(y,2))]);
subplot('Position',[0.07 0.41 0.4 0.19])
imagesc(y);
title(getString(message('signal:edr:OriginalSignal','Y')));
xlim([1 max(size(x,2),size(y,2))]);
subplot('Position',[0.07 0.11 0.4 0.19]);
m = size(x,2);
n = size(y,2);
% give equal weight to both images where they overlap
if m>n
 imagesc([(y+x(:,1:n))/2 x(:,n+1:end)]);
else
 imagesc([(x+y(:,1:m))/2 y(:,m+1:end)]);
end
title(getString(message('signal:edr:OverlaidOriginals')));

% Plot aligned signals in second column of plots
subplot('Position',[0.54 0.71 0.4 0.19]);
imagesc(x(:,ix));
title(getString(message('signal:edr:AlignedSignal','X')));
subplot('Position',[0.54 0.41 0.4 0.19]);
imagesc(y(:,iy));
title(getString(message('signal:edr:AlignedSignal','Y')));
subplot('Position',[0.54 0.11 0.4 0.19]);
imagesc(x(:,ix)+y(:,iy));
title(getString(message('signal:edr:OverlaidAlignment')));

% Add/modify master title by inserting an (invisible) axes.
hTitle = findall(gcf,'Tag','dtw_title');
if isempty(hTitle)
  hAxes = axes('Position',[0.09 0.09 .86 .86], ...
               'Visible','off', 'Tag','dtw_axes');
  hTitle = get(hAxes,'Title');
end

set(hTitle, ...
  'Visible','on', ...
  'Tag','dtw_title', ...
  'String', sprintf('%s: %i', ...
    getString(message('signal:edr:Edit')), dist));
  
allax=findobj(gcf,'type','axes','-not','tag','dtw_axes');
for iax=1:numel(allax)
  axes(allax(iax)); %#ok<LAXES>
end
