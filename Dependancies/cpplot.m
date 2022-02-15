function  cpplot(y, statistic, icp, residue)
%CPPLOT Helper function for plotting changepoints
%   This file is for internal use only and may be removed in a future
%   release

%   Copyright 2017 The MathWorks, Inc.

% create a new plot
m = size(y,1);
if m<=1
  hAxes = setupSinglePlot;
  cpplot1D(hAxes, y, statistic, icp);
elseif m==2
  [hAxes,hAxes2] = setupDualPlot;
  cpplot2D(hAxes, y, statistic, icp);
  cpplot1D(hAxes2, y, statistic, icp);
elseif m==3
  [hAxes,hAxes2] = setupDualPlot;
  cpplot3D(hAxes, y, statistic, icp);
  cpplot1D(hAxes2, y, statistic, icp);
else
  [hAxes,hAxes2] = setupDualPlot;
  cpplotND(hAxes, y, icp);
  cpplot1D(hAxes2, y, statistic, icp);
  linkaxes([hAxes,hAxes2],'x');
end

% title the (top) plot
numStr = getString(message('signal:findchangepts:NumberOfChangepoints',num2str(numel(icp))));
if any(strcmp(statistic,{'std','rms'}))
  errStr = getString(message('signal:findchangepts:TotalLogWeightedDispersion',num2str(residue)));
else
  errStr = getString(message('signal:findchangepts:TotalResidualError',num2str(residue)));
end
title({numStr, errStr});
hAxes.Parent.NextPlot = 'replacechildren';

%-------------------------------------------------------------------------
function hAxes = setupSinglePlot
hF = gcf;
if ~isempty(hF.Children) && strcmp(hF.NextPlot,'add')
  hF = figure;
else
  clf(hF);
end
hAxes = axes(hF);

%-------------------------------------------------------------------------
function [hAxes1,hAxes2] = setupDualPlot
hF = gcf;
if ~isempty(hF.Children) && strcmp(hF.NextPlot,'add')
  hF = figure;
else
  clf(hF);
end
hAxes2 = subplot(2,1,2,'Parent',hF);
hAxes1 = subplot(2,1,1,'Parent',hF);

%-------------------------------------------------------------------------
function  cpplot1D(hAxes, y, statistic, icp)
m = size(y,1);
n = size(y,2);
if n>1
  xlim(hAxes,[1 n]);
end

K = length(icp);
nseg = K+1;

plot(y.','Parent',hAxes);
colors = lines(size(y,1));

% plot trend for each region
if ~strcmp(statistic,'rms')
  for r=1:m
    xData = [1 icp-0.5; icp-0.5 n; nan(1, nseg)];
    yData = NaN(3,nseg);
    istart = [1 icp];
    istop = [icp-1 n];
    if strcmp(statistic,'mean') || strcmp(statistic,'std')
      for s=1:nseg
        ix = (istart(s):istop(s))';
        yData(1:2,s) = mean(y(r,ix));
      end
    elseif n>1 % linear
      for s=1:nseg
        ix = (istart(s):istop(s));
        b1 = 0.5*(s>1);
        b2 = 0.5*(s<nseg);
        yData(1:2,s) = polyval(polyfit(ix,y(r,ix),1),[ix(1)-b1 ix(end)+b2]);
      end
    end
    xData = reshape(xData, 3*nseg, 1);
    yData = reshape(yData, 3*nseg, 1);
    if m==1
      line(xData,yData,'Color',hAxes.ColorOrder(2,:),'Parent',hAxes);
    else
      line(xData,yData,'Color',colors(r,:),'LineWidth',1.2,'Parent',hAxes);
    end
  end
end

yLim = hAxes.YLim;
% plot vertical bars at changepoints
xData = reshape([icp-0.5; icp-0.5; NaN(size(icp))],3*K,1);
yData = repmat([yLim NaN].',K,1);
line(xData,yData,'Parent',hAxes,'Color',hAxes.ColorOrder(5,:));
hAxes.YLim = yLim;

%-------------------------------------------------------------------------
function cpplot2D(hAxes, y, statistic, icp)
K = length(icp);
nseg = K+1;
istart = [1; icp(:)];
istop = [icp(:)-1; size(y,2)];

if strcmp(statistic,'linear')
  % superimpose best fit segments over line
  line(hAxes,y(1,:),y(2,:),'Color',1-hAxes.Color);
  xData = NaN(3,nseg);
  yData = NaN(3,nseg);
  for s=1:nseg
    ix = (istart(s):istop(s));
    b1 = 0.5*(s>1);
    b2 = 0.5*(s<nseg);
    xData(1:2,s) = polyval(polyfit(ix,y(1,ix),1),[ix(1)-b1; ix(end)+b2]);
    yData(1:2,s) = polyval(polyfit(ix,y(2,ix),1),[ix(1)-b1; ix(end)+b2]);
  end
  line(hAxes,xData,yData,'Color',hAxes.ColorOrder(4,:),'Marker','.','LineWidth',1.2);
else % just plot colored for each changepoint region
  colors = lines(nseg);
  for s=1:nseg
    ix = (istart(s):istop(s));
    line(hAxes,y(1,ix),y(2,ix),'Color',colors(s,:));
  end
end
hAxes.XColor = hAxes.ColorOrder(1,:);
hAxes.YColor = hAxes.ColorOrder(2,:);

%-------------------------------------------------------------------------
function cpplot3D(hAxes, y, statistic, icp)
K = length(icp);
nseg = K+1;
istart = [1; icp(:)];
istop = [icp(:)-1; size(y,2)];

if strcmp(statistic,'linear')
  % superimpose best fit segments over line
  plot3(hAxes,y(1,:),y(2,:),y(3,:),'Color',1-hAxes.Color);
  xData = zeros(2,nseg);
  yData = zeros(2,nseg);
  zData = zeros(2,nseg);
  for s=1:nseg
    ix = (istart(s):istop(s));
    b1 = 0.5*(s>1);
    b2 = 0.5*(s<nseg);
    xData(1:2,s) = polyval(polyfit(ix,y(1,ix),1),[ix(1)-b1; ix(end)+b2]);
    yData(1:2,s) = polyval(polyfit(ix,y(2,ix),1),[ix(1)-b1; ix(end)+b2]);
    zData(1:2,s) = polyval(polyfit(ix,y(3,ix),1),[ix(1)-b1; ix(end)+b2]);
  end
  line(hAxes,xData,yData,zData,'Color',hAxes.ColorOrder(4,:),'Marker','.','LineWidth',1.2);
else % just plot colored segments within each changepoint region
  colors = lines(nseg);
  for s=1:nseg
    ix = (istart(s):istop(s));
    line(hAxes,y(1,ix),y(2,ix),y(3,ix),'Color',colors(s,:));
  end
end
hAxes.XColor = hAxes.ColorOrder(1,:);
hAxes.YColor = hAxes.ColorOrder(2,:);
hAxes.ZColor = hAxes.ColorOrder(3,:);

%-------------------------------------------------------------------------
function cpplotND(hAxes, y, icp)
imagesc(y,'Parent',hAxes)
colorbar;
set(hAxes,'YDir','normal');
xData = [icp; icp; nan(1,length(icp))];
yData = repmat([ylim NaN]',1,length(icp));
line(hAxes,xData,yData,'Color','w','LineWidth',2)
