function findsignalplot(signal,data,nsignal,ndata,istart,istop,metric,maxseg,annotate)
%FINDSIGNALPLOT Helper function for plotting segments
%   This file is for internal use only and may be removed in a future
%   release

%   Copyright 2016 The MathWorks, Inc.

hF = gcf;
if isempty(hF.Children) || ...
    numel(hF.Children) == 1 && isa(hF.Children,'matlab.graphics.axis.Axes')
  hF.Tag = 'findsignalplot';
  hF.Position = [hF.Position(1:2)-[770 550]+hF.Position(3:4) 770 550];
elseif strcmp(hF.Tag,'findsignalplot')
  clf(hF);
  hF.Tag = 'findsignalplot';
end

if isvector(data) && isvector(signal)
  plotVectorSegments(signal, data, nsignal, ndata, istart, istop, maxseg, annotate);
else
  plotMatrixSegments(signal, data, nsignal, ndata, istart, istop, metric, annotate);
end

if strcmp(hF.Tag,'findsignalplot')
  hF.NextPlot = 'replace';
end

%-------------------------------------------------------------------------
function plotVectorSegments(signal,data,nsignal,ndata,istart,istop,maxseg,annotate)
if strcmp(annotate,'data')
  plotVectorData(data, istart, istop, maxseg);
elseif strcmp(annotate,'signal')
  subplot('Position',[0.1 0.1 0.15 0.8]);
  plotVectorSignal(signal);
  title(getString(message('signal:findsignal:Signal')));
  subplot('Position',[0.35 0.1 .6 .8]);
  plotVectorData(data, istart, istop, maxseg);
  title(getString(message('signal:findsignal:MatchesFound',numel(istart))));
elseif strcmp(annotate,'all')
  subplot('Position',[0.1 0.6 0.15 0.35]);
  plotVectorSignal(signal);
  title(getString(message('signal:findsignal:ActualSignal')));
  subplot('Position',[0.1 0.1 0.15 0.35]);
  plotVectorSignal(nsignal);
  title(getString(message('signal:findsignal:NormalizedSignal')));
  subplot('Position',[.35 0.6 0.6 0.35]);
  hLineData = plotVectorData(data, istart, istop, maxseg);
  title(getString(message('signal:findsignal:ActualData')));
  subplot('Position',[.35 0.1 0.6 0.35]);
  hLineNData = plotVectorData(ndata, istart, istop, maxseg);
  title(getString(message('signal:findsignal:NormalizedData')));
  if isreal(data) && ~isempty(hLineData) && ~isempty(hLineNData)
    linkaxes([hLineData(1).Parent, hLineNData(1).Parent],'x');
  end
end

%-------------------------------------------------------------------------
function hLine = plotVectorData(data, istart, istop, maxseg)
iseg = fetchSegmentVector(istart,istop);
if isreal(data)
  hLine = plot(1:numel(data),data,'-',iseg,data(iseg),'.');
  % add legend if not complex
  if ~isempty(istart)
    if numel(istart)==1
      legend(getString(message('signal:findsignal:Data')), ...
             getString(message('signal:findsignal:Signal')));
    else
      legend(getString(message('signal:findsignal:Data')), ...
             getString(message('signal:findsignal:Signals')));
    end
  else
    legend(getString(message('signal:findsignal:Data')))
  end    
else
  hLine = plot3(1:numel(data),real(data),imag(data),'-', ...
                iseg,real(data(iseg)),imag(data(iseg)),'.');
  ylabel('real');
  zlabel('imag');
  view(85,10);
end

if maxseg~=1
  title(getString(message('signal:findsignal:MatchesFound',numel(istart))));
end

%-------------------------------------------------------------------------
function plotMatrixSegments(signal,data,nsignal,ndata,istart,istop,metric,annotate)
if strcmp(annotate,'data')
  plotMatrixData(data, istart, istop, metric);
  title(getString(message('signal:findsignal:MatchesFound',numel(istart))));
elseif strcmp(annotate,'signal')
  subplot('Position',[0.05 0.1 .25 .8]);
  plotMatrixSignal(signal, metric);
  title(getString(message('signal:findsignal:Signal')));
  subplot('Position',[0.35 0.1 .6 .8]);
  plotMatrixData(data, istart, istop, metric);
  title(getString(message('signal:findsignal:MatchesFound',numel(istart))));
elseif strcmp(annotate,'all')
  subplot('Position',[0.05 0.6 0.3 0.35]);
  plotMatrixSignal(signal, metric);
  title(getString(message('signal:findsignal:ActualSignal')));
  subplot('Position',[0.05 0.1 0.3 0.35]);
  plotMatrixSignal(nsignal, metric);
  title(getString(message('signal:findsignal:NormalizedSignal')));
  subplot('Position',[0.4 0.6 0.55 0.35]);
  hImageData = plotMatrixData(data, istart, istop, metric);
  title(getString(message('signal:findsignal:ActualData')));
  subplot('Position',[0.4 0.1 0.55 0.35]);
  hImageNData = plotMatrixData(ndata, istart, istop, metric);
  title(getString(message('signal:findsignal:NormalizedData')));
  if ~isempty(hImageData) && ~isempty(hImageNData)
    linkaxes([hImageData(1).Parent, hImageNData(1).Parent],'x');
  end
end

%-------------------------------------------------------------------------
function plotVectorSignal(signal)
lineColor = get(0,'DefaultAxesColorOrder');
if isreal(signal)
  plot(signal,'.-','Color',lineColor(2,:));
  if numel(signal)>1
    xlim([1 numel(signal)]);
  else
    xlim([0 2]);
  end
else
  plot3(1:numel(signal),real(signal),imag(signal),'.-', ...
        'Color',lineColor(2,:));
  ylabel('real');
  zlabel('imag');
  view(85,10);
end

%-------------------------------------------------------------------------
function plotMatrixSignal(signal, metric)
if ~isreal(signal)
  % plot real portion on top, imaginary on bottom
  signal = [real(signal); imag(signal)];
end

if strcmp(metric,'symmkl')
  signal = 10*log10(signal);
end

imagesc(signal);
colorbar;

%-------------------------------------------------------------------------
function hImage = plotMatrixData(data, istart, istop, metric)
if ~isreal(data)
  % plot real portion on top, imaginary on bottom
  data = [real(data); imag(data)];
end

if strcmp(metric,'symmkl')
  data = 10*log10(data);
end

% plot just the segments if too many are found
if numel(istart) > 10
  iseg = fetchSegmentVector(istart,istop);
  sig = NaN(size(data));
  sig(:,iseg) = data(:,iseg);
  hImage = imagesc(sig);
  colorbar;
else
  % otherwise plot and mask unmatched data with semi-transparent patches
  
  % plot data
  hImage = imagesc(data);
  title(getString(message('signal:findsignal:Data')));
  colorbar;

  % create patches to mask unmatched data in consecutive index order
  [istart, idx] = sort(istart);
  istop = istop(idx);
  lastx = 0;
  
  % compute y-ordinates of patch with black transparancy
  % traverse corners in clockwise order [x2,y2; x1,y2; x1,y1; x2,y1]
  y = hImage.Parent.YLim([2 2 1 1]);
  c = [0 0 0];
  
  for i=1:numel(istart)
    if istart(i)>lastx+1
      x = [lastx+0.5 istart(i)-0.5 istart(i)-0.5 lastx+0.5];
      patch(x, y, c, 'EdgeColor','none','FaceAlpha',0.625);
    end
    lastx = istop(i);
  end
  
  if lastx<size(data,2)
    x = [lastx+0.5 size(data,2) size(data,2) lastx+0.5];
    patch(x, y, c, 'EdgeColor','none','FaceAlpha',0.625);
  end
end

%-------------------------------------------------------------------------
function iseg = fetchSegmentVector(istart,istop)
iseg = zeros(1,sum(istop-istart+1));
niseg = 0;
for i=1:numel(istart)
  npoints = istop(i)-istart(i)+1;
  iseg(niseg+1:niseg+npoints) = istart(i):istop(i);
  niseg=niseg+npoints;
end
