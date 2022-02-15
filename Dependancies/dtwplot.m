function dtwplot(x, y, ix, iy, dist, metric)
%DTWPLOT Helper function for plotting dynamically time warped signals
%
%   For vector X and Y, DTWPLOT plots the original signals and aligned
%   signals, where the original signals are placed in a (smaller) plot.
%
%   For matrix X and Y, DTWPLOT plots X and Y as images (via IMAGESC).  
%   To facilitate visual alignment of features, DTWPLOT additionally
%   plots X and Y together with 50% transparency where they overlap.

%   Copyright 2015 The MathWorks, Inc.

if isvector(x) && isvector(y)
  plotVectorAlignment(x, y, ix, iy, dist, metric);
else
  plotMatrixAlignment(x, y, ix, iy, dist, metric);
end

%-------------------------------------------------------------------------
function plotVectorAlignment(x, y, ix, iy, dist, metric)
% plot original signals in the first plot
if isreal(x)
  hax1 = subplot('Position',[0.13 0.74 0.78 0.17]); 
  plot(1:numel(x),x,'.-', ...
       1:numel(y),y,'.-')
else
  hax1 = subplot('Position',[0.13 .8 .75 .15]);
  plot3(1:numel(x),real(x),imag(x),'.-', ...
        1:numel(y),real(y),imag(y),'.-');
  view(-1.9,19.8);
  ylabel('real');
  zlabel('imag');
end
title(getString(message('signal:dtw:OriginalSignals')));

% plot aligned signals in the second plot
if isreal(x)
  hax2 = subplot('Position',[0.13 0.11 0.78 0.49]);
  plot(1:numel(ix),x(ix),'.-', ...
       1:numel(iy),y(iy),'.-');
else
  hax2 = subplot('Position',[0.13 0.11 .75 .55]);
  plot3(1:numel(ix),real(x(ix)),imag(x(ix)),'.-', ...
        1:numel(iy),real(y(iy)),imag(y(iy)),'.-');
  view(-3.5,13.4);
  ylabel('real');
  zlabel('imag');
end
titleStr = sprintf('%s (%s: %f)', ...
   getString(message('signal:dtw:AlignedSignals')), ...
   getString(message(['signal:dtw:' upper(metric(1)) metric(2:end)])), ...
   dist);
title(titleStr);
hax1.Tag = 'OriginalSignals';
hax2.Tag = 'AlignedSignals';


%-------------------------------------------------------------------------
function plotMatrixAlignment(x, y, ix, iy, dist, metric)
if ~isreal(x) || ~isreal(y)
  % split each matrix into its real and imaginary portions
  x = [real(x); imag(x)]; y = [real(y); imag(y)];
end

% Plot original signals in first column of plots
hax1 = subplot('Position',[0.07 0.71 0.4 0.19]);
imagesc(x);
title(getString(message('signal:dtw:OriginalSignal','X')));
xlim([1 max(size(x,2),size(y,2))]);
hax2 = subplot('Position',[0.07 0.41 0.4 0.19]);
imagesc(y);
title(getString(message('signal:dtw:OriginalSignal','Y')));
xlim([1 max(size(x,2),size(y,2))]);
hax3 = subplot('Position',[0.07 0.11 0.4 0.19]);
m = size(x,2);
n = size(y,2);
% give equal weight to both images where they overlap
if m>n
 imagesc([(y+x(:,1:n))/2 x(:,n+1:end)]);
else
 imagesc([(x+y(:,1:m))/2 y(:,m+1:end)]);
end
title(getString(message('signal:dtw:OverlaidOriginals')));

% Plot aligned signals in second column of plots
hax4 = subplot('Position',[0.54 0.71 0.4 0.19]);
imagesc(x(:,ix));
title(getString(message('signal:dtw:AlignedSignal','X')));
hax5 = subplot('Position',[0.54 0.41 0.4 0.19]);
imagesc(y(:,iy));
title(getString(message('signal:dtw:AlignedSignal','Y')));
hax6 = subplot('Position',[0.54 0.11 0.4 0.19]);
imagesc(x(:,ix)+y(:,iy));
title(getString(message('signal:dtw:OverlaidAlignment')));

% Add/modify master title by inserting an (invisible) axes.
% prevent unwanted interaction with zoom by turning off handle visibility.
hTitle = findall(gcf,'Tag','dtw_title');
if isempty(hTitle)
  hAxes = axes('Position',[0.09 0.09 .86 .86], ...
               'Visible','off', 'Tag', 'dtw_axes');
  hTitle = get(hAxes,'Title');
end

set(hTitle, ...
  'Visible','on', ...
  'Tag','dtw_title', ...
  'String', sprintf('%s: %f', ...
    getString(message(['signal:dtw:' upper(metric(1)) metric(2:end)])), ...
    dist));

hax1.Tag = 'OriginalSignalX';
hax2.Tag = 'OriginalSignalY';
hax3.Tag = 'OverlaidOriginals';
hax4.Tag = 'AlignedSignalX';
hax5.Tag = 'AlignedSignalY';
hax6.Tag = 'OverlaidAlignment';

allax=findobj(gcf,'type','axes','-not','tag','dtw_axes');
for iax=1:numel(allax)
  axes(allax(iax)); %#ok<LAXES>
end
