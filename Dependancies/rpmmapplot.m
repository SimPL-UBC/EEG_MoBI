function rpmmapplot(map, freqVect, timeVect, rpmVect, res, type, scale, amp, rpmRef, timeRef)
%RPMMAPPLOT Plotter function for RPM maps.
%   This function is for internal use only. It may be removed in the future.
%
%   Inputs:
%   map      - order or frequency map matrix
%   freqVect - vector of frequencies or orders
%   timeVect - time vector
%   rpmVect  - vector of RPM values
%   res      - frequency or order resolution
%   type     - map type, can be 'Order' or 'Frequency'
%   scale    - scale can be linear or dB
%   amp      - can be 'rms', 'peak' or 'power'

% Copyright 2015-2019 The MathWorks, Inc.

% Create and render UI Framework
hFrame = createUIFramework(res,type);

% Create body of the GUI (e.g., axes, images, etc)
createGraphics(hFrame, map, freqVect, timeVect, rpmVect, scale, amp, rpmRef, timeRef);
end

%----------------------------------------------------------
function toggle2D3DView(~,~,hfig)

if nargin < 3
  hfig=gcbf;
end
ud = get(hfig,'UserData');
haxSpec = ud.hax(1);

% If currently in 3D, then change to 2D and vice versa
if ud.View3DStatus
  % Go to 2D view
  set(haxSpec,'View',[0 90],...
    'XTickLabel','');
  set(get(haxSpec,'XLabel'),'String','');
  set(haxSpec, 'Position',ud.haxSpecSizeLarge);
  set(get(haxSpec,'YLabel'),'Rotation',90);  
  maxZ = max(get(haxSpec,'ZLim'));
  set(ud.hspecX,'ZData',[maxZ, maxZ]); 
  set(ud.hspecY,'ZData',[maxZ, maxZ]);  
  set(ud.hCursors3D,'Visible','off');
  set(haxSpec,'YTickMode','auto');
  set(haxSpec,'XTickMode','auto');
  
else
  % Go to 3D view
  set(haxSpec,'View',[45,50],...
    'XTickLabelMode','auto');
  set(haxSpec, 'Position',ud.haxSpecSizeMedium);   
  set(get(haxSpec,'YLabel'),'Rotation',0);
  
  % Move cursors to the floor of the 3D plot
  minZ = min(get(haxSpec,'ZLim'));
  set(ud.hspecX,'ZData',[minZ, minZ]);
  set(ud.hspecY,'ZData',[minZ, minZ]);  
  set(ud.hCursors3D,'Visible','on');
end


ud.View3DStatus = ~ud.View3DStatus;
set(hfig,'UserData',ud);

updateAxesWithEngUnits(hfig);
update3DControlItems(hfig);

end

%----------------------------------------------------------
function rotate3DViewLeft(~,~,hfig)

if nargin < 3
  hfig=gcbf;
end
ud = get(hfig,'UserData');
if ~ud.View3DStatus 
  % Not in 3D mode so return
  return;
end

haxSpec = ud.hax(1);
newView = get(haxSpec,'View');
newView(1) = min(75, newView(1)+5);

if newView(1) == 20 
  rot = 48;
elseif newView(1) == 15
  rot = 55;
else
  rot = 0;
end

if newView(1) >= 25 && newView(1) <= 65
  set(haxSpec, 'Position',ud.haxSpecSizeMedium);
  set(get(haxSpec,'YLabel'),'Rotation',rot);
  set(haxSpec,'View',newView);
elseif newView(1) >= 15 && newView(1) <= 75
  %set(haxSpec, 'Position',ud.haxSpecSizeSmall);
  set(haxSpec,'View',newView);
  set(get(haxSpec,'YLabel'),'Rotation',rot);
end

updateAxesWithEngUnits(hfig);

end

%----------------------------------------------------------
function rotate3DViewRight(~,~,hfig)

if nargin < 3
  hfig=gcbf;
end
ud = get(hfig,'UserData');
if ~ud.View3DStatus
  % Not in 3D mode so return  
  return;
end

haxSpec = ud.hax(1);
newView = get(haxSpec,'View');
newView(1) = max(15, newView(1)-5);

if newView(1) == 20 
  rot = 48;
elseif newView(1) == 15
  rot = 55;
else
  rot = 0;
end

if newView(1) >= 25 && newView(1) <= 65
  set(haxSpec, 'Position',ud.haxSpecSizeMedium);     
  set(haxSpec,'View',newView);  
  set(get(haxSpec,'YLabel'),'Rotation',rot);
  
elseif newView(1) >= 15 && newView(1) <= 75
  %set(haxSpec, 'Position',ud.haxSpecSizeSmall);
  set(haxSpec,'View',newView);
  set(get(haxSpec,'YLabel'),'Rotation',rot);
end

updateAxesWithEngUnits(hfig);

end

%----------------------------------------------------------
function reset3DView(~,~,hfig)
if nargin < 3
  hfig=gcbf;
end
ud = get(hfig,'UserData');
if ~ud.View3DStatus
  % Not in 3D mode so return
  return;
end
haxSpec = ud.hax(1);
set(haxSpec,'View',[45,50]);
set(haxSpec, 'Position',ud.haxSpecSizeMedium);
set(get(haxSpec,'YLabel'),'Rotation',0);

updateAxesWithEngUnits(hfig);

end

%----------------------------------------------------------
function toggleZoom(ctrl,~,hfig)
if nargin < 3
  hfig=gcbf;
end
ud = get(hfig,'UserData');

if isprop(ctrl,'Checked')
  % Menu item was clicked - if not checked, then user wants to turn
  % clicked mode on, if checked user wants to turn clicked mode off
  if strcmp(ctrl.Checked,'off')
    % Turn mode on
    ud.ZoomOn = true;    
    ud.ZoomMode = getZoomModeFromTag(ctrl);
  else
    % Turn mode off
    ud.ZoomOn = false;    
    ud.ZoomMode = [];
  end
elseif strcmp(ctrl.State,'on')
  % Toolbar button was pressed
  ud.ZoomOn = true;
  ud.ZoomMode = getZoomModeFromTag(ctrl);
else
  % Toolbar button was de-pressed
  ud.ZoomOn = false;
  ud.ZoomMode = [];  
end
  
set(hfig,'UserData',ud);

updateZoomControlItems(hfig);

end

%----------------------------------------------------------
function zoomFull(hco,eventStruct,hfig)

if nargin < 3
  hfig=gcbf;
end

ud = get(hfig,'UserData');

focusTimeReset(hco, eventStruct, hfig); % reset focus window to full extent

% Get time range of full spectrogram and update
surfXdata = get(ud.hsurface,'XData');
newXlim = [min(surfXdata), max(surfXdata)];
set(ud.hax(1),'XLim',newXlim);

surfYdata = get(ud.hsurface,'YData');
newYlim = [min(surfYdata), max(surfYdata)];
set(ud.hax(1),'YLim',newYlim);

% Update the thumbnail:
set(ud.hthumb,'XData',newXlim([1 2 2 1 1]));

updateAxesWithEngUnits(hfig);

end

%---------------------------------------------------------------
function setCrosshairs(hfig,x,y)

% Get current point in axis ASAP:
ud = get(hfig,'UserData');
hax = ud.hax(1);

% Update cache
if isempty(x)
  x = ud.crosshair.xctr;
else
  ud.crosshair.xctr = x;
end

if nargin == 3
  if isempty(y)
    y = ud.crosshair.yctr;
  else
    ud.crosshair.yctr = y;
  end
else
  y = ud.crosshair.yctr;
end
set(hfig,'UserData',ud);

% Update crosshairs
set([ud.hspecY ud.htimeY], 'XData',[x x]);
set(ud.hspecX, 'YData',[y y]);

% Update 3D cursors
zLim = get(hax, 'ZLim');
z0 = zLim(1);
z1 = zLim(2);

xData = get(ud.hspecX,'XData');
x0 = xData(1);
x1 = xData(2);
xPatch = [x0 x1 x1 x0];
yPatch = [y   y  y  y];
zPatch = [z0 z0 z1 z1];
set(ud.hCursors3D(1), 'XData',xPatch, 'YData',yPatch, 'ZData', zPatch);

yData = get(ud.hspecY,'YData');
y0 = yData(1);
y1 = yData(2);

xPatch = [x   x  x  x];
yPatch = [y0 y1 y1 y0];
zPatch = [z0 z0 z1 z1];
set(ud.hCursors3D(2), 'XData',xPatch, 'YData',yPatch, 'ZData', zPatch);

% Update readouts
updateTimeReadout(hfig);
updateFreqReadout(hfig);
updateAmpReadout(hfig);
updateRpmReadout(hfig);

end

%---------------------------------------------------------------
function centerCursors(~,~,hfig)

if nargin < 3
  hfig=gcbf;
end

ud = get(hfig,'UserData');

% Determine center of spectrogram axis:
xlim=get(ud.hax(1),'XLim');
ylim=get(ud.hax(1),'YLim');

setCrosshairs(hfig,mean(xlim),mean(ylim));

end

%---------------------------------------------------------------
function bringCursorsIntoZoomRegion(hfig)
% If cursors are already in zoom region do nothing, but if they are not,
% then bring them to center of region.

if nargin < 3
  hfig=gcbf;
end

ud = get(hfig,'UserData');

% Determine center of spectrogram axis:
xlim=get(ud.hax(1),'XLim');
ylim=get(ud.hax(1),'YLim');

xPos = [];
yPos = [];
if ud.crosshair.yctr < ylim(1) || ud.crosshair.yctr > ylim(2) 
  yPos = mean(ylim);
end
if ud.crosshair.xctr < xlim(1) || ud.crosshair.xctr > xlim(2) 
  xPos = mean(xlim);
end

setCrosshairs(hfig,xPos,yPos);

end


%---------------------------------------------------------------
function wbmotionThumb(~,~)
% thumbnail motion

% Get current point in axis ASAP:
hfig   = gcbf;
ud     = get(hfig,'UserData');
hax    = ud.hax(2);
cp     = get(hax,'CurrentPoint');
currX  = cp(1,1);

xmotion = currX - ud.thumb.origPt;
width   = ud.thumb.width;
xdata   = ud.thumb.xdata + xmotion;

% Constrain to axis limits, so we don't lose cursor:
xlim=get(hax,'XLim');
minXdata = min(xdata);
maxXdata = max(xdata);
if minXdata < xlim(1)
  xdata=[xlim(1) xlim([1 1])+width xlim([1 1])];
elseif maxXdata > xlim(2)
  xdata = [xlim(2)-width xlim([2 2]) xlim([2 2])-width];
end

% If the patch is larger than the zoom window
if min(xdata)<=xlim(1) && max(xdata)>=xlim(2)
  return
end

% Update the thumbnail:
set(ud.hthumb,'XData',xdata);

% Scroll the spectrogram:
set(ud.hax(1),'XLim',xdata(1:2));  % [1 3]

if ud.View3DStatus
  % Update spectrogram time axis if in 3d view
  updateSpecTimeAxesWithEngUnits(hfig);
end

end

%---------------------------------------------------------------
function wbmotionThumbLeft(~,~)
% thumbnail LEFT motion

% Get current point in axis ASAP:
hfig   = gcbf;
ud     = get(hfig,'UserData');

% current object may be either the patch, or the signal line
hax    = ud.hax(2);
cp     = get(hax,'CurrentPoint');
currX  = cp(1,1);

xmotion = currX - ud.thumb.origPt;
xdata   = ud.thumb.xdata;
xdata([1 4 5]) = xdata([1 4 5]) + xmotion;

% Constrain to axis limits, so we don't lose cursor:
xlim = get(hax,'XLim');
minXdata = min(xdata);
if minXdata < xlim(1)
  xdata=[xlim(1) xdata([2 3])' xlim([1 1])];
elseif minXdata >= xdata(2)
  xdata = ud.thumb.xdata;
end

% Update the thumbnail:
set(ud.hthumb,'XData',xdata);

% Scroll the spectrogram:
set(ud.hax(1),'XLim',xdata(1:2));

if ud.View3DStatus
  % Update spectrogram time axis if in 3d view
  updateSpecTimeAxesWithEngUnits(hfig);
end

end

%---------------------------------------------------------------
function wbmotionThumbright(~,~)
% thumbnail RIGHT motion

% Get current point in axis ASAP:
hfig   = gcbf;
ud     = get(hfig,'UserData');

hax    = ud.hax(2);
cp     = get(hax,'CurrentPoint');
currX = cp(1,1);

xmotion = currX - ud.thumb.origPt;
xdata   = ud.thumb.xdata;
xdata([2 3]) = xdata([2 3]) + xmotion;

% Constrain to axis limits, so we don't lose cursor:
xlim = get(hax,'XLim');
maxXdata = max(xdata);
if maxXdata > xlim(2)
  xdata(2:3) = xlim(2);
elseif maxXdata <= xdata(1)
  xdata = ud.thumb.xdata;
end

% Update the thumbnail:
set(ud.hthumb,'XData',xdata);

% Scroll the spectrogram:
set(ud.hax(1),'XLim',xdata(1:2));

if ud.View3DStatus
  % Update spectrogram time axis if in 3d view
  updateSpecTimeAxesWithEngUnits(hfig);
end

end

%---------------------------------------------------------------
function wbupThumb(~,~)

% set spectrogram and time-slice xlims
hfig = gcbf;
ud = get(hfig,'UserData');

% Commented out, due to flash:
%     set(ud.hax([1 3]),'XLim',xlim);
%
% This is fine:
%     set(ud.hax(1),'XLim',xlim);
%
% the following line causes flash in the spect image
% this is the time-slice axis:
% why does this affect the spectrogram axis? (overlap? clipping?)
%     set(ud.hax(3),'XLim',xlim);

changePtr(gcbf,'hand');
installCursorFcns(gcbf,'thumb');

% Turn back on image axis visibility,
% which was turned off during wbdownThumb
% so that it does not flash while panning:
% Leave off!
%set(ud.hax(1),'Visible','on');
%
% Turn on crosshair visibility, which was shut off
% during thumbnail panning (wbdownThumb):
set([ud.hspecY ud.hspecX],'Visible','on');

if ud.View3DStatus
  % Update spectrogram time axis if in 3d view
  updateSpecTimeAxesWithEngUnits(hfig);
end

end

%---------------------------------------------------------------
function wbupThumbLeft(~,~)
% set spectrogram and time-slice xlims
hfig = gcbf;
ud = get(hfig,'UserData');
xdata = get(ud.hthumb,'XData');
xlim = [min(xdata) max(xdata)];

% Commented out, due to flash:
%set(ud.hax([1 3]),'XLim',xlim);
%
% This is fine:
set(ud.hax(1),'XLim',xlim);
%
% the following line causes flash in the spect image
% this is the time-slice axis:
%set(ud.hax(3),'XLim',xlim);

changePtr(gcbf,'ldrag');
installCursorFcns(gcbf,'thumbleft');

% Turn on crosshair visibility, which was shut off during thumbnail panning
% (wbdownThumb):
set([ud.hspecY ud.hspecX],'Visible','on');

if ud.View3DStatus
  % Update spectrogram time axis if in 3d view
  updateSpecTimeAxesWithEngUnits(hfig);
end

end

%---------------------------------------------------------------
function wbupThumbRight(~,~)
% set spectrogram and time-slice xlims
hfig = gcbf;
ud = get(hfig,'UserData');
xdata = get(ud.hthumb,'XData');
xlim = [min(xdata) max(xdata)];

% Commented out, due to flash:
%set(ud.hax([1 3]),'XLim',xlim);
%
% This is fine:
set(ud.hax(1),'XLim',xlim);
%
% the following line causes flash in the spect image
% this is the time-slice axis:
%set(ud.hax(3),'XLim',xlim);

changePtr(gcbf,'rdrag');
installCursorFcns(gcbf,'thumbright');

% Turn on crosshair visibility, which was shut off during thumbnail panning
% (wbdownThumb):
set([ud.hspecY ud.hspecX],'Visible','on');

if ud.View3DStatus
  % Update spectrogram time axis if in 3d view
  updateSpecTimeAxesWithEngUnits(hfig);
end

end

%---------------------------------------------------------------
function updateStatus(hfig, str)
% updateStatus Update status text.

ud = get(hfig,'UserData');
hstatusbar = ud.hStatusBar;
if(~isempty(hstatusbar))
  if nargin == 2 && ischar(str)
    hstatusbar.Text = str;
  elseif nargin == 1
    strres = [ud.ResTxt, ud.ResValue ' ' ud.ResUnits];    
    ud.hStatus.Text = strres;
  else
    return;
  end
end
end

%---------------------------------------------------------------
function setCmapLimits(hfig, new_dr)
% Set new colormap limits

ud = get(hfig,'UserData');
haxSpec = ud.hax(1);
haxCbar = ud.hax(3);
himageCbar = ud.himageCbar;

% Set new dynamic range limits into spectrogram image
set(haxSpec,'CLim',new_dr);

% colorbar is 1:256
% Actual spectrogram dynamic range is orig_dr new spectrogram dynamic range
% is new_dr
origDr = get(himageCbar,'YData');
diffDr = new_dr - origDr;
cmapIndices_per_dB = 256./diff(origDr);  % a constant
diffClim = diffDr .* cmapIndices_per_dB;
cbarClim = [1 256] + diffClim;
set(himageCbar,'CDataMapping','scaled');  % do during creation
set(haxCbar,'CLim',cbarClim);

end

%---------------------------------------------------------------
function resetCmapLimits(~,~)
% Reset colormap limits to dynamic range of spectrogram data

hfig = gcbf;
ud = get(hfig,'UserData');
origDr = get(ud.himageCbar,'YData');
setCmapLimits(hfig, origDr);

end

%---------------------------------------------------------------
function manualCmapLimits(~,~,hfig)
% manualCmapLimits Manual change to colormap dynamic range limits

if nargin < 3
  hfig = gcbf;
end
ud = get(hfig,'UserData');
haxSpec = ud.hax(1);

% Prompt for changes to cmap limits:
clim = get(haxSpec,'CLim');
  
prompt={getStrFromCat('topColor'),getStrFromCat('bottomColor')};

def = {num2str(clim(2)), num2str(clim(1))};
dlgTitle = getStrFromCat('color');
lineNo=1;
strs=inputdlg(prompt,dlgTitle,lineNo,def);
if isempty(strs)
  return
end
new_dr = [str2double(strs{2}) str2double(strs{1})];

setCmapLimits(hfig,new_dr);

end

%---------------------------------------------------------------
function wbmotionCmap(~,~)
% wbmotionCmap Graphical change to colormap dynamic range limits

hfig = gcbf;
ud = get(hfig,'UserData');
haxSpec = ud.hax(1);
haxCbar = ud.hax(3);

% Determine cursor starting and current points ASAP:
cp    = get(haxCbar,'CurrentPoint');
newPt = cp(1,2);  % y-coord only
dy    = newPt - ud.cbar.origPt;

% if SelectionType was normal, update top or bottom of colorbar, only,
% depending on whether user started drag in the top or bottom of bar,
% respectively.
%
% if SelectionType was extend, update both top AND bottom of bar
% simultaneously, translating colormap region.
if strcmp(ud.cbar.SelectionType,'extend')
  change_dr = [dy dy];
else
  if ud.cbar.StartInTop
    change_dr = [0 dy];
  else
    change_dr = [dy 0];
  end
end
newDr = ud.cbar.starting_dr + change_dr;
if diff(newDr)<=0
  newDr = ud.cbar.starting_dr;
end

% Colorbar range is 1 to 256. Actual spectrogram dynamic range is orig_dr.
% New spectrogram dynamic range is new_dr.
origDr = get(ud.himageCbar,'YData');    % a constant
cmapIndices_per_dB = 256./diff(origDr);  % a constant
diffDr = newDr - origDr;
diffClim = diffDr .* cmapIndices_per_dB;
cbarClim = [1 256] + diffClim;

if diff(cbarClim)>0
  % Protect against poor choice of values
  set(haxCbar,'CLim',cbarClim,'UserData',newDr);
end

% We defer setting the new dynamic range limits into the spectrogram image
% axis, as it will create too much flash. Instead, on button-up, the new
% limit is set.  See wbupCmap() for details.

% Set new dynamic range limits into spectrogram image. Note: userdata could
% be empty if this is the first entry...
set(haxSpec,'CLim',newDr);

end

%---------------------------------------------------------------
function isChange = changePtr(hfig, newPtr)

% Get current pointer name:
ud = get(hfig,'UserData');

% Is this a change in pointer type?
isChange = ~strcmp(ud.currPtr,newPtr);
if isChange
  setptr(hfig, newPtr);
  ud.currPtr = newPtr;
  set(hfig,'UserData',ud);
end

end

%---------------------------------------------------------------
function wbmotionGeneral(~,~,hfig)
% General button motion
%
% Determines if cursor is over a crosshair. If so, changes pointer and
% installs crosshair buttondowns If not, changes back to normal cursor and
% general buttondowns as necessary.

if nargin < 3
  hfig = gcbf;
end
ud = get(hfig,'UserData');

[isOverHV, isSpecgramAxis, isCmapAxis,isTopHalfCmap, isThumb] = overCrosshair(hfig);

if ~any(isOverHV)
  % Not hovering over a crosshair
  
  if isSpecgramAxis
    % Do nothing over the spectrogram axes unless in zoom XY mode
    if ud.ZoomOn
      if strcmp(ud.ZoomMode,'ZoomOut')
        changePtr(hfig, 'glassminus');
        installCursorFcns(hfig,'zoomout');
      else
        changePtr(hfig, 'glassplus');
        installCursorFcns(hfig,'zoomdrag');
      end
      
    else
      changePtr(hfig,'arrow');
      installCursorFcns(hfig,'general');
    end
    
  elseif isCmapAxis
    % Over the colormap axes
    % Install the up/down pointer:
    if isTopHalfCmap
      if changePtr(hfig,'udrag')
        updateStatus(hfig,getStrFromCat('upperShiftToTranslate'));
        installCursorFcns(hfig,'cmap');
      end
    else
      if changePtr(hfig,'ddrag')
        updateStatus(hfig,getStrFromCat('lowerShiftToTranslate'));
        installCursorFcns(hfig,'cmap');
      end
    end
    
  elseif any(isThumb)
    % Over thumbnail - isThumb is a 3-element vector, [left center right],
    %  indicating whether cursor is over left edge, right edge, or is over
    %  the general thumbnail patch itself.
    
    % Install appropriate pointer:
    if isThumb(1)
      % Over left edge
      if changePtr(hfig,'ldrag')
        installCursorFcns(hfig,'thumbleft');
      end
    elseif isThumb(3)
      % Over right edge
      if changePtr(hfig,'rdrag')
        installCursorFcns(hfig,'thumbright');
      end
    else
      % Over general patch region
      if changePtr(hfig,'hand')
        installCursorFcns(hfig,'thumb');
      end
    end
    
  else
    % Not over a special axes:
    if changePtr(hfig,'arrow')
      installCursorFcns(hfig,'general');
    end
  end
  
else
  % Pointer is over a crosshair (vert or horiz or both)
  if all(isOverHV)
    % Over both horiz and vert (near crosshair center):
    if changePtr(hfig,'fleur')
      installCursorFcns(hfig,'hvcross');
    end
  elseif isOverHV(1)
    % Over H crosshair
    if changePtr(hfig,'uddrag')
      installCursorFcns(hfig,'hcross');
    end
  else
    % Over V crosshair
    if changePtr(hfig,'lrdrag')
      installCursorFcns(hfig,'vcross');
    end
  end
end

end

%---------------------------------------------------------------
function [y,isSpecAxis,isCmapAxis,isTopHalfCmap, isThumb] = overCrosshair(hfig)
% Is the cursor hovering over the crosshairs?
%
% There are two crosshairs, one an H-crosshair, the other a V-crosshair.
% The H and V crosshairs span several different axes.
%
% Function returns a 2-element vector, indicating whether the cursor is
% currently over the H- and/or V-crosshairs.
%
% y = [isOverH isOverV]

y             = [0 0];
isSpecAxis = 0;
isCmapAxis    = 0;
isTopHalfCmap = 0;
isThumb       = [0 0 0];  % left, middle, right regions

% First, are we over any axes?
hax = overAxes(hfig);
if isempty(hax)
  return;
end  % not over an axis

% Get current point in axis:
cp = get(hax,'CurrentPoint');
ud = get(hfig,'UserData');

isCmapAxis    = (hax == ud.hax(3));
isSpecAxis = (hax == ud.hax(1));

% Determine if any horiz or vert crosshairs are in this axis - store as
% [anyHoriz anyVert]:
hasHVCrossHairs = [hax == ud.hax(1) any(hax == ud.hax(1:2))];

% Is cursor in colormap axis?
if (isCmapAxis)
  % is cursor in top half of colormap axis?
  orig_dr = get(ud.hax(1),'CLim');
  isTopHalfCmap = (cp(1,2) >= sum(orig_dr)/2);
end

if any(hasHVCrossHairs)
  % Get cursor & crosshair positions:
  crosshairPos = [ud.crosshair.xctr ud.crosshair.yctr];
  curretPt     = cp(2,1:2); 
  cursorDelta  = abs(crosshairPos - curretPt);
  axisDx       = diff(get(hax,'XLim'));
  axisDy       = diff(get(hax,'YLim'));
  axisDelta    = [axisDx axisDy];
  
  % Is cursor within 1 percent of crosshair centers?
  % Limit test uses the reciprocal of the percentage tolerance
  %   1-percent -> 1 / 0.01 = 100
  % 1.5-percent -> 1 / 0.015 ~= 67
  %   2-percent -> 1 / 0.02 = 50
  %
  % Finally, allow a true result only if the axis
  % has a crosshair of the corresponding type
  %
  y = fliplr(cursorDelta * 67 < axisDelta) & hasHVCrossHairs;
end

% Are we over the thumbnail patch?
% Check if we're over the time axis:
if (hax == ud.hax(2))
  % Get thumb patch limits:
  xdata=get(ud.hthumb,'XData');
  xlim=[min(xdata) max(xdata)];
  
  % Is cursor over general patch area?
  thumb_delta = xlim - cp(1,1);
  isThumb(2) = thumb_delta(1)<=0 & thumb_delta(2)>=0;
  
  % Is cursor over left or right thumbnail edge?
  % Use same tolerance as crosshair test:
  axisDx        = diff(get(hax,'XLim'));
  isThumb([1 3]) = (abs(thumb_delta) * 67 < axisDx);
end

end 

%---------------------------------------------------------------
function h = overAxes(hfig)
% overAxes Determine if pointer is currently over an axis of the figure;
% the axis list comes from the figure UserData (ud.hax).

fig = get(hfig);
ud = fig.UserData;

obj = hittest(hfig);

switch obj.Tag
  case {'rpmmapplotSpecAxesTag','rpmmapplotSurfaceTag',...
      'rpmmapplotTopOfSpecAxesTag','rpmmapplotSpecCursorXTag',...
      'rpmmapplotSpecCursorYTag','rpmmapplotSpecCursorPatch1Tag',...
      'rpmmapplotSpecCursorPatch2Tag'}
    h = ud.hax(1);
    
  case {'rpmmapplotPannerAxesTag','rpmmapplotPannerPatchTag',...
          'rpmmapplotPannerLineTag','rpmmapplotPannerMarkerTag',...
          'rpmmapplotPannerCursorTag'}
    
    h = ud.hax(2);
    
  case {'rpmmapplotCBarAxesTag','rpmmapplotCBarImgTag'}
    h = ud.hax(3);
        
  otherwise
    h = [];
end

end
%---------------------------------------------------------------
function y = isLeftClick(hfig)

% Keywords for key/button combinations:
%         Left    Right
%   none: normal  alt
%  Shift: extend  alt
%   Ctrl: alt     alt
% Double: open    alt

y=strcmp(get(hfig,'SelectionType'),'normal');
end

%---------------------------------------------------------------
function wbdownHCross(~,~)
% window button down in h-crosshair mode
if ~isLeftClick(gcbf)
  return;
end
ud = get(gcbf,'UserData');
set(ud.hCursors3D(1),'FaceAlpha',.4);
installCursorFcns(gcbf,'hcross_buttondown');
wbmotionCross([],[],'h');
end

%---------------------------------------------------------------
function wbdownVCross(~,~)
% window button down in v-crosshair mode
if ~isLeftClick(gcbf)
  return;
end
ud = get(gcbf,'UserData');
set(ud.hCursors3D(2),'FaceAlpha',.4);
installCursorFcns(gcbf,'vcross_buttondown');
wbmotionCross([],[],'v');
end

%---------------------------------------------------------------
function wbdownHVCross(~,~)
% window button down in hv-crosshair mode
if ~isLeftClick(gcbf)
  return;
end
ud = get(gcbf,'UserData');
set(ud.hCursors3D,'FaceAlpha',.4);
installCursorFcns(gcbf,'hvcross_buttondown');
wbmotionCross([],[],'hv');

end

%---------------------------------------------------------------
function wbdownThumb(~,~)
% window button down in thumbnail mode
if ~isLeftClick(gcbf)
  return;
end

% cache y-coord of pointer
ud = get(gcbf,'UserData');
haxTime = ud.hax(2);
cp = get(haxTime,'CurrentPoint');
xdata = get(ud.hthumb,'XData');
width = max(xdata)-min(xdata);

ud.thumb.origPt = cp(1,1);   % x-coord only
ud.thumb.width  = width;
ud.thumb.xdata  = xdata;
set(gcbf,'UserData',ud);

changePtr(gcbf,'closedhand');
installCursorFcns(gcbf,'thumb_buttondown');

% Turn off image axis visibility,
% so that it does not flash while panning:
%
% off permanently now:
%set(ud.hax(1),'Visible','off');
%
% Turn off crosshair visibility:
set([ud.hspecY ud.hspecX],'Visible','off');

end

%---------------------------------------------------------------
function wbdownThumbLeft(~,~)

% window button down in LEFT thumbnail mode
if ~isLeftClick(gcbf)
  return;
end

% cache y-coord of pointer
ud = get(gcbf,'UserData');
hax_time = ud.hax(2);
cp = get(hax_time,'CurrentPoint');
xdata = get(ud.hthumb,'XData');
width = max(xdata)-min(xdata);

ud.thumb.origPt = cp(1,1);   % x-coord only
ud.thumb.width  = width;
ud.thumb.xdata  = xdata;
set(gcbf,'UserData',ud);

installCursorFcns(gcbf,'thumbleft_buttondown');

% Turn off crosshair visibility:
set([ud.hspecY ud.hspecX],'Visible','off');

end

%---------------------------------------------------------------
function wbdownThumbRight(~,~)

% window button down in LEFT thumbnail mode
if ~isLeftClick(gcbf)
  return;
end

% cache y-coord of pointer
ud = get(gcbf,'UserData');
haxTime = ud.hax(2);
cp = get(haxTime,'CurrentPoint');
xdata = get(ud.hthumb,'XData');
width = max(xdata)-min(xdata);

ud.thumb.origPt = cp(1,1);   % x-coord only
ud.thumb.width  = width;
ud.thumb.xdata  = xdata;
set(gcbf,'UserData',ud);

installCursorFcns(gcbf,'thumbright_buttondown');

% Turn off crosshair visibility:
set([ud.hspecY ud.hspecX],'Visible','off');

end

%----------------------------------------------------
function wbdownCmap(~,~)
% window button down in colormap mode

hfig = gcbf;

% Only allow left (normal) or shift+left (extend)
st = get(hfig,'SelectionType');
i = find(strncmp(st,{'normal','extend','open'}, length(st)));
if isempty(i)
  return;
end

if i==3
  % open dynamic range menu
  manualCmapLimits([],[],hfig);
  return
elseif i==2
  % Shift+left button = translate,
  % show up/down cursor during drag
  % NOTE: cannot update cursor when shift is pressed
  %       but no mouse button is pressed (no event!)
  changePtr(hfig,'uddrag');
end

ud = get(hfig,'UserData');

% cache y-coord of pointer
haxCbar = ud.hax(3);
cp = get(haxCbar,'CurrentPoint');
ud.cbar.origPt = cp(1,2);   % y-coord only
ud.cbar.SelectionType = st; % normal or extend

% The current clim is in the spectrogram image. We want to know the midpoint
% of this
origDr = get(ud.hax(1),'CLim');
ud.cbar.midPt = sum(origDr)/2;

% Determine if pointer went down in top or bottom half of colorbar:
ud.cbar.StartInTop = (ud.cbar.origPt >= ud.cbar.midPt);

% Cache original dynamic range:
haxSpec = ud.hax(1);
ud.cbar.starting_dr = get(haxSpec,'CLim');
set(hfig,'UserData',ud);

installCursorFcns(hfig,'cmap_buttondown');

% Set initial clim into userdata in case motion callback not performed
% (motion updates userdata). wbupCmap reads the userdata
%
% Turn off visibility during drag to prevent flash
set(haxCbar, ...
  'UserData',ud.cbar.starting_dr, ...
  'Visible','off');
end

%---------------------------------------------------------------
function wbupHCross(~,~)
% window button up in h-crosshair mode
installCursorFcns(gcbf,'hcross');
ud = get(gcbf,'UserData');
set(ud.hCursors3D,'FaceAlpha',.05);
end

%---------------------------------------------------------------
function wbupVCross(~,~)
% window button up in v-crosshair mode
installCursorFcns(gcbf,'vcross');
ud = get(gcbf,'UserData');
set(ud.hCursors3D,'FaceAlpha',.05);
end

%---------------------------------------------------------------
function wbupHVCross(~,~)
% window button up in hv-crosshair mode
installCursorFcns(gcbf,'hvcross');
ud = get(gcbf,'UserData');
set(ud.hCursors3D,'FaceAlpha',.05);
end

%---------------------------------------------------------------
function wbupCmap(~,~)
% window button up in colormap mode
installCursorFcns(gcbf,'cmap');

% Set new dynamic range limits into spectrogram image. Note: userdata could
% be empty if this is the first entry...
ud = get(gcbf,'UserData');
haxCbar = ud.hax(3);
set(ud.hax(1),'CLim',get(haxCbar,'UserData'));
set(haxCbar,'Visible','on'); % re-enable axis vis

% Set new status msg, since it doesn't update in the installCursorFcns fcn
% for cmap callbacks Do this by calling the general mouse-motion fcn:
wbmotionGeneral([],[]);

end
%---------------------------------------------------------------
function wbdownZoom(~,~)

if ~isLeftClick(gcbf)
  return;
end
hfig = gcbf;
ud = get(hfig,'UserData');
haxSpec = ud.hax(1);
cp = get(haxSpec,'CurrentPoint');
ud.ZoomOrigPt = cp(2,1:2);
ud.ZoomDestPt = [];

if ud.View3DStatus
  set(ud.hthumbZoom3D,'Xdata',NaN,'YData',NaN,'ZData',NaN,'Visible','on');
else
  set(ud.hthumbZoom,'Xdata',NaN,'YData',NaN,'Visible','on');
end
installCursorFcns(gcbf,'zoom_buttondown');

% Turn off crosshair visibility:
set([ud.hspecY ud.hspecX],'Visible','off');
set(ud.hCursors3D,'Visible','off');

set(hfig,'UserData',ud);
end

%---------------------------------------------------------------
function wbmotionZoom(~,~,hfig)
% Paint thumb region based on mouse movement and cache zoom limits
if nargin < 3
  hfig = gcbf;
end
ud = get(hfig,'UserData');
haxSpec = ud.hax(1);
axXlims = haxSpec.XLim;
axYlims = haxSpec.YLim;
cp = get(haxSpec,'CurrentPoint');

xo = ud.ZoomOrigPt(1);
yo = ud.ZoomOrigPt(2);

if cp(1,2) < 0    
  % Do not move the region if pointer is outside axes
  return;
end

switch ud.ZoomMode
  case 'ZoomIn'
    x1 = cp(2,1);
    y1 = cp(2,2);
  case 'ZoomX'
    x1 = cp(2,1);
    yo = min(axYlims);
    y1 = max(axYlims);    
  case 'ZoomY'
    y1 = cp(2,2);
    xo = min(axXlims);
    x1 = max(axXlims);        
  case 'ZoomOut'
    return;
end
  
ud.ZoomDestPt = [x1,y1];

minX = max(min(xo,x1),min(axXlims));
maxX = min(max(xo,x1),max(axXlims));

minY = max(min(yo,y1),min(axYlims));
maxY = min(max(yo,y1),max(axYlims));

ud.Zoom_XLimits = [minX maxX];
ud.Zoom_YLimits = [minY maxY];

xCoords = [minX maxX maxX minX];
yCoords = [minY minY maxY maxY];

if ud.View3DStatus
  % Draw 3D zoom region
  x0 = minX;
  x1 = maxX;
  y0 = minY;
  y1 = maxY;
  
  zLim = get(haxSpec, 'ZLim');
  z0 = zLim(1);
  z1 = zLim(2);
  
  x = [x0 x1 x1 x0; x1 x1 x1 x0; x1 x1 x0 x0; x0 x1 x0 x0];
  y = [y0 y0 y1 y1; y0 y0 y1 y0; y0 y1 y1 y0; y0 y1 y1 y1];
  z = [z0 z0 z0 z0; z0 z1 z1 z0; z1 z1 z1 z1; z1 z0 z0 z1];
  
  set(ud.hthumbZoom3D,'Xdata',x,'YData',y,'ZData',z);  
else
  % Draw 2D zoom patch
  maxZ = max(get(haxSpec,'ZLim'));
  set(ud.hthumbZoom,'Xdata',xCoords,'YData',yCoords,'ZData',[maxZ maxZ maxZ maxZ]);
end

set(hfig,'UserData',ud);

end
%---------------------------------------------------------------
function wbupZoom(~,~)

hfig = gcbf;
ud = get(hfig,'UserData');

% Remove zoom thumb regions
set(ud.hthumbZoom,'Xdata',NaN,'YData',NaN,'Visible','off');
set(ud.hthumbZoom3D,'Xdata',NaN,'YData',NaN,'ZData',NaN,'Visible','off');

% Turn on crosshair visibility:
set([ud.hspecY ud.hspecX],'Visible','on');
set(ud.hCursors3D,'Visible','on');

switch ud.ZoomMode
  case 'ZoomIn'
    % Go back to zoom mode
    installCursorFcns(gcbf,'zoomdrag');
    doZoomX(hfig,ud);
    doZoomY(hfig,ud);
  case 'ZoomX'
    installCursorFcns(gcbf,'zoomdrag');
    doZoomX(hfig,ud)  
  case 'ZoomY'
    installCursorFcns(gcbf,'zoomdrag');
    doZoomY(hfig,ud)            
  case 'ZoomOut'
    installCursorFcns(gcbf,'zoomout');
    doZoomX(hfig,ud,true);
    doZoomY(hfig,ud,true);  
end

% Update ticks and labels
updateAxesWithEngUnits(hfig);
bringCursorsIntoZoomRegion(hfig);

end

%---------------------------------------------------------------
function doZoomX(hfig,ud,zoomOutFlag)

if nargin == 1
  ud = get(hfig,'UserData');
end
if nargin < 3
  zoomOutFlag = false;
end

% Pixel spacing of data. x axis points can be non uniformly sampled so pick
% smallest width (values at end of data).
surfXdata = get(ud.hsurface,'XData');
minXData = min(surfXdata);
maxXData = max(surfXdata);

dataDx = surfXdata(end)-surfXdata(end-1);
currentXLimits = get(ud.hax(1),'XLim');

if zoomOutFlag
  % Do 50% zoom out in x
  newDeltaXLim = diff(currentXLimits)*2;
  deltaXData = maxXData - minXData;
  if newDeltaXLim > deltaXData
    % New delta limits exceed data limits so do full zoom
    zoomXLim = [minXData, maxXData];
  else
    xCenter = ud.ZoomOrigPt(1);
    xLeft = xCenter - newDeltaXLim/2;
    xRight = xCenter + newDeltaXLim/2;
    if xLeft < minXData
      % Min limits are smaller than min data so pick left half of data
      % limits
      zoomXLim = [minXData, minXData+deltaXData];
    elseif xRight > maxXData
      % Max limits are larger than max data so pick right half of data
      % limits
      zoomXLim = [maxXData-deltaXData, maxXData];
    else
      zoomXLim = [xLeft, xRight];
    end
  end

elseif isempty(ud.ZoomDestPt) || diff(ud.Zoom_XLimits) == 0
  % Single click so do 50% zoom in x
  
  % If we can place the center at the point where mouse was clicked without
  % exceeding min/max data limits do so. Otherwise if mouse was clicked to
  % the left of the current axes center, zoom to the left half or if mouse
  % was clicked to the right of the current axes center, zoom to the left
  % half.   
  newDeltaXLim = diff(currentXLimits)/2;
  xCenter = ud.ZoomOrigPt(1);
  xLeft = xCenter - newDeltaXLim/2;
  xRight = xCenter + newDeltaXLim/2;
  if xLeft >minXData && xRight < maxXData
    zoomXLim = [xLeft, xRight];    
  else
    % If left of center zoom to left half, otherwise, zoom to right half
    xCenter = mean(currentXLimits);
    if ud.ZoomOrigPt(1) < xCenter
      zoomXLim = [currentXLimits(1), xCenter];
    else
      zoomXLim = [xCenter, currentXLimits(2)];
    end    
  end
  if diff(zoomXLim) <= 0
    % Do not zoom if zoom delta is zero
    return;
  end

else
  % Check that the zoom region is not zero width in X. If it is, then
  % increase the zoom width to half pixel width.
  
  % Get X zoom limits
  zoomXLim = ud.Zoom_XLimits;
  
  % Zoom delta X
  dxLim = diff(zoomXLim);
  
  if dxLim <= 0
    % Increase second zoom X limit by 0.5*pixel width or max current limit
    zoomXLim(2) = min(currentXLimits(2),zoomXLim(2)+dataDx/2);
  end     
end

% Perform zoom (also update the panner thumbnail x data). Do this only if
% we have different min max values,
if diff(zoomXLim) > 0
  set(ud.hax(1),'XLim',zoomXLim);
  set(ud.hthumb,'XData',zoomXLim([1 2 2 1 1]));
end

end
%---------------------------------------------------------------
function doZoomY(hfig,ud,zoomOutFlag)

if nargin == 1
  ud = get(hfig,'UserData');
end
if nargin < 3
  zoomOutFlag = false;
end

% Pixel spacing of data.
surfYdata = get(ud.hsurface,'YData');
minYData = min(surfYdata);
maxYData = max(surfYdata);
dataDy = surfYdata(2)-surfYdata(1);
currentYLimits = get(ud.hax(1),'YLim');

if zoomOutFlag
  % Do 50% zoom out in y
  newDeltaYLim = diff(currentYLimits)*2;
  deltaYData = maxYData - minYData;
  if newDeltaYLim > deltaYData
    % New delta limits exceed data limits so do full zoom
    zoomYLim = [minYData, maxYData];
  else
    yCenter = ud.ZoomOrigPt(2);
    yLow = yCenter - newDeltaYLim/2;
    yUp = yCenter + newDeltaYLim/2;    
    if yLow < minYData
      % Min limits are smaller than min data so pick lower half of data
      % limits
      zoomYLim = [minYData, minYData+deltaYData];      
    elseif yUp > maxYData
      % Max limits are larger than max data so pick upper half of data
      % limits      
      zoomYLim = [maxYData-deltaYData, maxYData];   
    else
      zoomYLim = [yLow, yUp];
    end
  end

elseif isempty(ud.ZoomDestPt) || diff(ud.Zoom_YLimits) == 0
  % Single click so do 50% zoom in x
  
  % If we can place the center at the point where mouse was clicked without
  % exceeding min/max data limits do so. Otherwise if mouse was clicked
  % above the current axes center, zoom to the upper half or if mouse was
  % clicked below the current axes center, zoom to the lower half.  
  newDeltaYLim = diff(currentYLimits)/2;
  yCenter = ud.ZoomOrigPt(2);
  yLow = yCenter - newDeltaYLim/2;
  yUp = yCenter + newDeltaYLim/2;
  if yLow > minYData && yUp < maxYData
    zoomYLim = [yLow, yUp];
  else    
    % If below of center zoom to bottom half, otherwise, zoom to top half
    yCenter = mean(currentYLimits);
    if ud.ZoomOrigPt(2) < yCenter
      zoomYLim = [currentYLimits(1), yCenter];
    else
      zoomYLim = [yCenter, currentYLimits(2)];
    end
  end
  if diff(zoomYLim) <= 0
    % Do not zoom if zoom delta is zero
    return;
  end
else
  % Check that the zoom region is not zero width in Y. If it is, then
  % increase the zoom width to half pixel width.
  
  % Get Y zoom limits
  zoomYLim = ud.Zoom_YLimits;
  
  % Zoom delta Y
  dyLim = diff(zoomYLim);
  
  if dyLim <= 0
    % Increase second zoom Y limit by 0.5*pixel width or max current limit
    zoomYLim(2) = min(currentYLimits(2),zoomYLim(2)+dataDy/2);
  end    
end

% Perform zoom only if we have 2 different values
if diff(zoomYLim) > 0
  set(ud.hax(1),'YLim',zoomYLim);
end

end

%---------------------------------------------------------------
function [i,j] = getAdjustedCrosshairIdx(hfig)
% Find image matrix coordinate pair (j,i) under crosshair.
% Adjust crosshair for "half-pixel offset" implicit in image display

ud = get(hfig,'UserData');
xc = ud.crosshair.xctr;
yc = ud.crosshair.yctr;
hsurface = ud.hsurface;
im = get(hsurface,'CData');

% Get image pixel size:
xdata = get(hsurface,'XData');
if length(xdata)>1
  dx = xdata(2)-xdata(1);
else
  dx=0;
end

ydata=get(hsurface,'YData');
if length(ydata)>1
  dy = ydata(2)-ydata(1);
else
  dy=0;
end

% Remove half a pixel size from apparent cursor Position:
xc = xc-dx/2;
yc = yc-dy/2;

% Find pixel coordinate under the crosshair:
i = find(xc>=xdata);
if isempty(i)
  i=1;
else
  i=i(end)+1;
end
j=find(yc>=ydata);
if isempty(j)
  j=1;
else
  j=j(end)+1;
end
sz=size(im);
if i>sz(2)
  i=sz(2);
end
if j>sz(1)
  j=sz(1);
end

end

%---------------------------------------------------------------
function v = getSpecVal(hfig)

ud         = get(hfig,'UserData');
surfCData  = get(ud.hsurface,'CData');
[i,j]      = getAdjustedCrosshairIdx(hfig);
v          = double(surfCData(j,i));

end

%---------------------------------------------------------------
function v = getRPMVal(hfig)

ud    = get(hfig,'UserData');
yd    = get(ud.htimePlot,'YData');
i     = getAdjustedCrosshairIdx(hfig);
v     = double(yd(i));

end

%---------------------------------------------------------------
function updateTimeReadout(hfig)

ud = get(hfig,'UserData');
t = ud.crosshair.xctr;
prefix='T = ';

% Update time readout
[y,~,u] = engunits(t, 'latex','time');
if strcmp(u, 'secs')
  u = 's';
end

str = [prefix sprintf('%.3f',y) ' ' u];
set(ud.htextTime,'String',str);

end

%---------------------------------------------------------------
function updateFreqReadout(hfig)

ud = get(hfig,'UserData');
isOrderMap = strcmpi(ud.MapType, 'Order');
f = ud.crosshair.yctr;

if isOrderMap    
  str= [getStrFromCat('orderReadoutLbl') ' = '  sprintf('%.3f',f)];
else
  [y,~,u] = engunits(f,'latex');
  prefix='F = ';
  unitsPostFix = 'Hz';
  str=[prefix sprintf('%.3f',y) ' ' u unitsPostFix];
end

% Update freq readout
set(ud.htextFreq,'String',str);

end

%---------------------------------------------------------------
function updateRpmReadout(hfig)

ud = get(hfig,'UserData');
prefix = 'RPM = ';

rpm = getRPMVal(hfig);
str = [prefix sprintf('%.3f',rpm)];
set(ud.htextRpm,'String',str);

end

%---------------------------------------------------------------
function updateAmpReadout(hfig)

ud = get(hfig,'UserData');

switch ud.AmplitudeType
  case 'rms'
    prefix = [getStrFromCat('rmsAmpReadoutLbl') ' = '];
  case 'peak'
    prefix = [getStrFromCat('peakAmpReadoutLbl') ' = '];
  case 'power'
    prefix = [getStrFromCat('power') ' = '];
end

if strcmpi(ud.Scale, 'dB')
  postfix = ' dB';
else
  postfix = '';
end

% Update amplitude readout
a = getSpecVal(hfig);
str=[prefix sprintf('%.3f',a) postfix];
set(ud.htextAmp,'String',str);

end

%---------------------------------------------------------------
function wbmotionCross(~,~,sel)
% motion callback during horiz/vert-crosshair selection
% sel='h', 'v', or 'hv'

% Get current point in axis ASAP:
hfig = gcbf;
hco  = gco;
switch get(hco,'Type')
  case 'axes'
    hax=hco;
  otherwise
    hax=get(hco,'Parent');
end

cp    = get(hax,'CurrentPoint');

if cp(1,2) < 0    
  % Do not move the cursors if pointer is outside axes
  return;
end

ud    = get(hfig,'UserData');
x     = cp(2,1);
y     = cp(2,2);

switch sel
  case 'h'
    x = ud.crosshair.xctr;
  case 'v'
    y = ud.crosshair.yctr;
end 

% Constrain to axis limits, so we don't lose cursor:
if any(sel == 'v')
  xlim = get(hax,'XLim');
  if x < xlim(1)
    x = xlim(1);
  elseif x > xlim(2)
    x = xlim(2);
  end
end

if any(sel == 'h')
  ylim = get(hax,'YLim');
  if y < ylim(1)
    y = ylim(1);
  elseif y > ylim(2)
    y = ylim(2);
  end
end

setCrosshairs(hfig,x,y);
end

%---------------------------------------------------------------
function installCursorFcns(hfig,cursorType)

switch lower(cursorType)
  case 'none'
    dn     = [];
    motion = [];
    up     = [];
    status = '';
    
  case 'general'
    dn     = [];
    motion = @wbmotionGeneral;
    up     = [];
    status = '';
    
  case 'thumb'
    % button not pushed, thumbnail highlighted
    dn     = @wbdownThumb;
    motion = @wbmotionGeneral;
    up     = [];
    status = getStrFromCat('panZoomWindow');
    
  case 'thumb_buttondown'
    % button pushed, thumbnail highlighted
    dn     = [];
    motion = @wbmotionThumb;
    up     = @wbupThumb;
    status = getStrFromCat('releaseToSetZoomWindow');
    
  case 'thumbleft'
    % button not pushed, left thumbnail edge highlighted
    dn     = @wbdownThumbLeft;
    motion = @wbmotionGeneral;
    up     = [];
    status = getStrFromCat('adjustZoomWinLeft');
    
  case 'thumbleft_buttondown'
    % button pushed, thumbnail highlighted
    dn     = [];
    motion = @wbmotionThumbLeft;
    up     = @wbupThumbLeft;
    status = getStrFromCat('releaseToSetZoomWindow');
    
  case 'thumbright'
    % button not pushed, right thumbnail edge highlighted
    dn     = @wbdownThumbRight;
    motion = @wbmotionGeneral;
    up     = [];
    status = getStrFromCat('adjustZoomWinRight');
    
  case 'thumbright_buttondown'
    % button pushed, right thumbnail edge highlighted
    dn     = [];
    motion = @wbmotionThumbright;
    up     = @wbupThumbRight;
    status = getStrFromCat('releaseToSetZoomWindow');
    
  case 'hcross'
    % button not pushed, h-crosshair highlighted
    dn     = @wbdownHCross;
    motion = @wbmotionGeneral;
    up     = [];
    status = getStrFromCat('moveHCursor');
    
  case 'hcross_buttondown'
    % button pushed while over horiz cross-hair
    dn     = [];
    motion = {@wbmotionCross,'h'};
    up     = @wbupHCross;
    status = getStrFromCat('releaseToUpdateCursor');
    
  case 'vcross'
    dn     = @wbdownVCross;
    motion = @wbmotionGeneral;
    up     = [];
    status = getStrFromCat('moveVCursor');
    
  case 'vcross_buttondown'
    dn     = [];
    motion = {@wbmotionCross,'v'};
    up     = @wbupVCross;
    status = getStrFromCat('releaseToUpdateCursor');
    
  case 'hvcross'
    dn     = @wbdownHVCross;
    motion = @wbmotionGeneral;
    up     = [];
    status = getStrFromCat('moveCrossHairCursor');
    
  case 'hvcross_buttondown'
    dn     = [];
    motion = {@wbmotionCross,'hv'};
    up     = @wbupHVCross;
    status = getStrFromCat('releaseToUpdateCursor');
    
    % Change dynamic range of colormap
  case 'cmap'
    dn     = @wbdownCmap;
    motion = @wbmotionGeneral;
    up     = [];
    % Status is set in wbmotionGeneral function since it depends on which
    % pointer we're using
    status = -1;
    
  case 'cmap_buttondown'
    dn     = [];
    motion = @wbmotionCmap;
    up     = @wbupCmap;
    status =  getStrFromCat('releaseToUpdateCmap');
    
  case 'zoomdrag'
    dn     = @wbdownZoom;
    motion = @wbmotionGeneral;
    up     = [];
    status =  getStrFromCat('zoomDragStatus');
    
  case 'zoomout'
    dn     = @wbdownZoom;
    motion = @wbmotionGeneral;
    up     = [];
    status =  getStrFromCat('zoomOutStatus');
    
  case 'zoom_buttondown'
    dn     = [];
    motion = @wbmotionZoom;
    up     = @wbupZoom;
    status = getStrFromCat('zoomButtonUpStatus');
    
  otherwise
    error(message('signal:rpmmapplot:invalidParamCursorFcn'));
end

set(hfig, ...
  'WindowButtonDownFcn',  dn, ...
  'WindowButtonMotionFcn',motion, ...
  'WindowButtonUpFcn',    up)

updateStatus(hfig,status);

end

%---------------------------------------------------------------
function resizeFig(~,~)
% Callback to resize the figure

updateAxesWithEngUnits(gcbf);

end

%---------------------------------------------------------------
function updateAxesWithEngUnits(hfig)
% Update the tick marks for axes that are using engineering units For
% example, a resize could have added or removed ticks, and the axes would
% no longer have the proper tick marks.

ud = get(hfig,'UserData');
%hFrame = getappdata(hfig, 'UIMgr');
if strcmp(hfig.Visible, 'on')
  haxSpec = ud.hax(1);
  haxTime = ud.hax(2);
  
  % Update freq-axis labels for engineering units
  yy = get(haxSpec,'YTick');
  [cs,eu] = convert2engstrs(yy);
  set(haxSpec,'YTickLabel',cs);
      
  if strcmpi(ud.MapType, 'Order')
    m = getMultiplier(eu);
    if isempty(m)
      set(get(haxSpec,'YLabel'),'String',getStrFromCat('orderReadoutLbl'));
    else
      set(get(haxSpec,'YLabel'),'String',...
        [getStrFromCat('orderReadoutLbl') ' ('  m ')']);
    end
  else   
    set(get(haxSpec,'YLabel'),...
      'String',[getStrFromCat('frequencyLbl') ' (' eu 'Hz)']);
  end
  
  if ud.View3DStatus
    % Set time axis label and tick labels for spectrogram plot
    yy = get(haxSpec,'XTick');
    [cs,eu] = convert2engstrs(yy,'time');
    set(haxSpec,'XTickLabel',cs);
    if strcmp(eu,'secs')
      eu = 's';
    end
    set(get(haxSpec,'XLabel'),'String',[getStrFromCat('timeLbl') ' (' eu ')']);
  end
    
  % Update time-axis labels for engineering units
  xt = get(haxTime,'XTick');
  [cs,eu] = convert2engstrs(xt,'time');
  set(haxTime,'XTickLabel',cs);
  if strcmp(eu,'secs')
    eu = 's';
  end
  set(get(haxTime,'XLabel'),'String',[getStrFromCat('timeLbl') ' (' eu ')']);
  
  if ud.RPMVect(end) >= 10e3
    yt = get(haxTime,'YTick');
    [cs,eu] = convert2engstrs(yt);
    set(haxTime,'YTickLabel',cs);
    m = getMultiplier(eu);
    if isempty(m)
      set(get(haxTime,'YLabel'),'String','RPM');
    else
      set(get(haxTime,'YLabel'),'String',...
        ['RPM' ' ('  m ')']);
    end
  end
end
end
%---------------------------------------------------------------
function updateSpecTimeAxesWithEngUnits(hfig)

ud = get(hfig,'UserData');

haxSpec = ud.hax(1);

% Set time axis label and tick labels for spectrogram plot
yy = get(haxSpec,'XTick');
[cs,eu] = convert2engstrs(yy,'time');
set(haxSpec,'XTickLabel',cs);
if strcmp(eu,'secs')
  eu = 's';
end
set(get(haxSpec,'XLabel'),'String',[getStrFromCat('timeLbl') ' (' eu ')']);

end

%---------------------------------------------------------------
function updateGUI(~, ~, hfig)

if nargin < 3
  hfig=gcbf;
end

ptr.ptr = get(hfig,'Pointer');
ptr.shape = get(hfig,'PointerShapeCData');
ptr.hot = get(hfig,'PointerShapeHotSpot');
setptr(hfig,'watch');  % set user's expectations...

ud = get(hfig,'UserData');
haxSpec     = ud.hax(1);
haxTime     = ud.hax(2);
haxCbar     = ud.hax(3);

b = ud.Map;
t = ud.t;
f = ud.f;

% Delete map matrix to reduce memory usage
ud.Map = [];

if strcmpi(ud.Scale, 'dB')
  % Handle -inf's:
  i_inf = find(isinf(b(:)));
  if ~isempty(i_inf)
    % Set all -inf points to next-lowest value:
    b(i_inf)=inf;
    min_val=min(b(:));
    b(i_inf)=min_val;
  end
end

blim = [min(b(:)) max(b(:))];
rpmlim = [min(ud.RPMRef) max(ud.RPMRef)];
% Leave 20% space in upper limit
deltaLim = diff(blim);
rpmdeltaLim = diff(rpmlim);
% If the upper and lower limits are the same, set rpmdeltaLim to the 
% rpm limit value 
if rpmdeltaLim == 0
  rpmdeltaLim = rpmlim(1); 
end

specZlim = [blim(1) blim(2)+deltaLim*.2];
specXlim = [0 max(t)];
specYlim = [0 max(f)];
% Leave 10% space in upper limit and lower limit
rpmYlim = [rpmlim(1)-rpmdeltaLim*0.1 rpmlim(2)+rpmdeltaLim*0.1];

% Update spectrogram
set(ud.hsurface,'CData',b, 'ZData',b, 'XData',t, 'YData',f);
haxSpec.View = [0 90];
set(haxSpec,'XLim',specXlim, 'YLim', specYlim, 'ZLim', specZlim);

% Update colorbar
set(ud.himageCbar,'XData',[0 1], 'YData',blim, 'CData', (1:256)');
set(haxCbar,'XLim',[0 1],'YLim',blim);

% This creates the RPM trace that goes on the pannable time plot below the
% main plot
set(ud.hrpmRef,'XData',ud.TimeRef,'YData',ud.RPMRef);
set(ud.htimePlot,'XData',t,'YData',ud.RPMVect);
set(haxTime,'XLim',specXlim);
set(haxTime,'YLim',rpmYlim);

% setup thumbnail patch
axylim = get(haxTime,'YLim');
ymax = axylim(2);
ymin = axylim(1);
tmax = max(t);
tmin = min(t);
set(ud.hthumb, ...
  'XData',[tmin tmax tmax tmin tmin], ...
  'YData',[ymin ymin ymax ymax ymin]);

% Reset crosshair positions
crosshair      = ud.crosshair;
crosshair.xctr = mean(specXlim);
crosshair.yctr = mean(specYlim);
timeYlim       = get(haxTime,'YLim');

% Crosshairs:
set(ud.hspecX, ...
  'XData',specXlim, ...
  'YData',[crosshair.yctr crosshair.yctr],...
  'ZData',[specZlim(2), specZlim(2)]);
set(ud.hspecY, ...
  'XData',[crosshair.xctr crosshair.xctr],...
  'YData',specYlim,...
  'ZData',[specZlim(2), specZlim(2)]);
set(ud.htimeY, ...
  'XData',[crosshair.xctr crosshair.xctr],...
  'YData',timeYlim);  

% Update user data:
ud.crosshair = crosshair;
set(hfig,'UserData',ud);

% Text readouts:
updateTimeReadout(hfig);
updateFreqReadout(hfig);
updateAmpReadout(hfig);
updateRpmReadout(hfig);

% Re-establish pointer cursor, etc:
set(hfig,'Pointer',ptr.ptr, ...
  'PointerShapeCData',ptr.shape, ...
  'PointerShapeHotSpot',ptr.hot);

%Set to normal display
zoomFull(0, 0, hfig);

end

%---------------------------------------------------------------
function printdlgCb(~,~)
printdlg(gcbf);
end

%---------------------------------------------------------------
function printpreviewCb(~,~, ~)
printpreview(gcbf);
end

%---------------------------------------------------------------
function closeCb(~,~)
delete(gcbf);
end

%---------------------------------------------------------------
function createGraphics(hfig, map, freqVect, timeVect, rpmVect, scale, amp, rpmRef, timeRef)
%createGraphics Render the graphics.

hfig.Colormap = parula;

ud = get(hfig,'UserData');

hVisParent = hfig;
if isappdata(hfig, 'UIMgr')
    hVisParent = get(getappdata(hfig, 'UIMgr'), 'hVisParent');
else
    hVisParent = hfig.findobj('Tag','VisualizationAreaContainer');
end

haxSpecSizeLarge  = [.1 .32 .74 .63];
haxSpecSizeMedium  = [.1 .37 .67 .58];
haxSpecSizeSmall = [.115 .38 .64 .57];

haxCbarSizeLarge = [.85 .321 .03 .630];
haxTimeSizeLarge = [.1 .15 .74 .14];

% specgram axes
haxSpec = axes('Parent',hVisParent, 'Position',haxSpecSizeLarge);
set(haxSpec, 'Box','on', 'XTickLabel','', 'XGrid', 'on','YGrid', 'on', ...
  'ZGrid','on','Tag','rpmmapplotSpecAxesTag');

if strcmpi(ud.MapType, 'order')
  str = getStrFromCat('orderMap');
else
  str = getStrFromCat('freqMap');
end
set(get(haxSpec,'Title'),'String', str);

% specgram surface
hSurf = surface('Parent',haxSpec, 'edgecolor', 'none',...
  'Tag','rpmmapplotSurfaceTag');

% colorbar
cmapLen = 256;
haxCbar = axes('Parent',hVisParent,'Position',haxCbarSizeLarge);
himageCbar = image([0 1],[0 1],(1:cmapLen)','Tag','rpmmapplotCBarImgTag');
set(himageCbar,'CDataMapping','scaled');
set(haxCbar, ...
  'Box','on', ...  
  'YDir','normal', ...
  'XTickLabel','', ...
  'XTick',[],...
  'YAxisLocation','right',...
  'Tag','rpmmapplotCBarAxesTag');

str = getStrFromCat(amp);
if strcmpi(scale,'dB')
  str = [str ' (dB)'];
end
set(get(haxCbar,'YLabel'),'String', str);

% RPM panner
haxTime = axes('Parent',hVisParent, 'Position',haxTimeSizeLarge);
% RPM plot - shows original RPM profile curve
htimePlot = line('Parent', haxTime, 'Color','b','Tag','rpmmapplotPannerMarkerTag');
set(htimePlot,'LineStyle','none','Marker','o','MarkerSize',2,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
% RPM markers - show RPM values at center of each computed spectrogram
% window
hrpmRef = line('Parent', haxTime, 'Color','b','Tag','rpmmapplotPannerLineTag');

set(haxTime, ...
  'Box','on',...
  'YAxisLocation','left',...
  'XGrid', 'on',...
  'YGrid', 'on',...
  'Tag','rpmmapplotPannerAxesTag');
axis tight
ylabel('Parent',haxTime, 'RPM'); % xlabel will be updated later

% thumbnail patch for panner
bgclr = 'b';
hthumb = patch([0 0 1 1 0], [0 1 1 0 0], bgclr, 'Parent',haxTime,...
  'Tag','rpmmapplotPannerPatchTag');
transparency = 0.1; % [0,1] controls transparency of patch FaceColor.
bgclr = [.1 .1 .1];  % FaceColor of the patch
set(hthumb,'EdgeColor',bgclr,'EdgeAlpha',transparency, ...
  'FaceColor',bgclr,'FaceAlpha',transparency);

% Put all spectrogram cursors, patches and zoom patches on an invisible
% axes on top of spectrogram axes. Link View, Position, scale and limits of
% spectrogram and cursor axes.
haxTopOfSpec = axes('Parent',hVisParent, 'Position',haxSpecSizeLarge);
set(haxTopOfSpec, 'Visible','off','Tag','rpmmapplotTopOfSpecAxesTag');

% haxSpec is first on vector so that its property values are set initially
% on linked objects. 
hLinkAxesProps = linkprop([haxSpec, haxTopOfSpec],{'Position','XLim',...
  'YLim','ZLim','XScale','YScale','View'});

% Cursors crosshair lines
hspecX = line('Parent',haxTopOfSpec,'Tag','rpmmapplotSpecCursorXTag');
hspecY = line('Parent',haxTopOfSpec,'Tag','rpmmapplotSpecCursorYTag');
htimeY = line('Parent',haxTime,'Tag','rpmmapplotPannerCursorTag');

% Cursors cross-patches for 3D plot
transparency = 0.05; % [0,1] controls transparency of patch FaceColor.
bgclr = [.5 .5 .5];  % FaceColor of the patch
hCursors3D = [...
  patch(NaN, NaN, bgclr, 'Parent',haxTopOfSpec,'Visible','off','Tag','rpmmapplotSpecCursorPatch1Tag'),...
  patch(NaN, NaN, bgclr, 'Parent',haxTopOfSpec,'Visible','off','Tag','rpmmapplotSpecCursorPatch2Tag')];
set(hCursors3D,'EdgeColor',bgclr,'EdgeAlpha',transparency, ...
  'FaceColor',bgclr,'FaceAlpha',transparency,'Visible','off');

% Thumbnail patch for zoom
hthumbZoom = patch([0 0 1 1 0], [0 1 1 0 0], bgclr, 'Parent',haxTopOfSpec);
transparency = 0.3; % [0,1] controls transparency of patch FaceColor.
bgclr = [.5 .5 .5];  % FaceColor of the patch
set(hthumbZoom,'EdgeColor',bgclr,'EdgeAlpha',transparency, ...
  'FaceColor',bgclr,'FaceAlpha',transparency,'Visible','off');

% Patches for  for 3D zoom
transparency = 0.2; % [0,1] controls transparency of patch FaceColor.
bgclr = [.5 .5 .5];  % FaceColor of the patch
hthumbZoom3D = [...
  patch([0 0 1 1 0], [0 1 1 0 0], bgclr, 'Parent',haxTopOfSpec,'Visible','off'),...
  patch([0 0 1 1 0], [0 1 1 0 0], bgclr, 'Parent',haxTopOfSpec,'Visible','off'),...
  patch([0 0 1 1 0], [0 1 1 0 0], bgclr, 'Parent',haxTopOfSpec,'Visible','off'),...
  patch([0 0 1 1 0], [0 1 1 0 0], bgclr, 'Parent',haxTopOfSpec,'Visible','off')];
set(hthumbZoom3D,'EdgeColor',bgclr,'EdgeAlpha',transparency, ...
  'FaceColor',bgclr,'FaceAlpha',transparency,'Visible','off');

% Text readouts:
w = 0.247;
h = 0.04;
xpos = 0.006;
ypos = 0.007;
lblXpos = .018;
lblYpos = 0.5;
font = {'FontSize',8};
haxReadout1 = axes('Parent',hVisParent,'Position',[xpos ypos w h],'Visible','off',...
  'Tag','rpmmapplotReadout1Tag');
patch([0 1 1 0 0],[0 0 1 1 0],'w');
htextFreq = text('Parent',haxReadout1, 'Position',[lblXpos lblYpos], font{:});

haxReadout2 = axes('Parent',hVisParent,'Position',[xpos+w ypos w h],'Visible','off',...
  'Tag','rpmmapplotReadout2Tag');
patch([0 1 1 0 0],[0 0 1 1 0],'w');
htextRpm = text('Parent',haxReadout2, 'Position',[lblXpos lblYpos], font{:});

haxReadout3 = axes('Parent',hVisParent,'Position',[xpos+2*w ypos w h],'Visible','off',...
  'Tag','rpmmapplotReadout3Tag');
patch([0 1 1 0 0],[0 0 1 1 0],'w');
htextTime = text('Parent',haxReadout3, 'Position',[lblXpos lblYpos], font{:});

haxReadout4 = axes('Parent',hVisParent,'Position',[xpos+3*w ypos w h],'Visible','off',...
  'Tag','rpmmapplotReadout4Tag');
patch([0 1 1 0 0],[0 0 1 1 0],'w');
htextAmp = text('Parent',haxReadout4, 'Position',[lblXpos lblYpos], font{:});

set([haxTopOfSpec haxCbar haxTime], 'SortMethod','childorder');

% Retain info in figure userdata
ud.hfig         = hfig;
ud.hax          = [haxSpec, haxTime, haxCbar];
ud.hspecX       = hspecX;
ud.hspecY       = hspecY;
ud.htimeY       = htimeY;
ud.htimePlot    = htimePlot;
ud.hrpmRef      = hrpmRef;
ud.htextTime    = htextTime;
ud.htextFreq    = htextFreq;
ud.htextAmp     = htextAmp;
ud.htextRpm     = htextRpm;
ud.htextStatus  = [];
ud.crosshair    = [];
ud.hsurface     = hSurf;
ud.himageCbar   = himageCbar;
ud.hthumb       = hthumb;
ud.hthumbZoom   = hthumbZoom;
ud.hthumbZoom3D = hthumbZoom3D;
ud.hCursors3D   = hCursors3D;
ud.currPtr      = '';

% Need to save link object or properties will not be linked
ud.hLinkAxesProps = hLinkAxesProps; 

ud.haxSpecSizeLarge  = haxSpecSizeLarge;
ud.haxSpecSizeMedium = haxSpecSizeMedium;
ud.haxSpecSizeSmall  = haxSpecSizeSmall;

ud.haxCbarSizeLarge  = haxCbarSizeLarge;
ud.haxTimeSizeLarge  = haxTimeSizeLarge;

ud.f             = freqVect;
ud.t             = timeVect;
ud.Map           = map;
ud.RPMVect       = rpmVect;
ud.RPMRef        = rpmRef;
ud.TimeRef       = timeRef;
ud.Scale         = scale;
ud.AmplitudeType = amp;

ud.ZoomOn = false;
ud.ZoomMode = [];
ud.ZoomOrigPt = [];
ud.ZoomDestPt = [];
ud.Zoom_XLimits = [];
ud.Zoom_YLimits = [];

ud.param_dlg = [];

set(hfig,'UserData',ud);

% Protect GUI from user plots, etc:
set([hfig ud.hax],'HandleVisibility','Callback');

% After GUI has all elements in it, install context help:
installContextHelp(hfig);
installContextMenus(hfig);

% Populate GUI with data, limits, etc:
updateGUI([],[],hfig);

% Enable general (non-segmenting) mouse functions:
installCursorFcns(hfig,'general');
set(hfig,'Visible','on');
if isappdata(hfig, 'UIMgr')
  %hFrame = getappdata(hfig, 'UIMgr');
  %hFrame.Enable = 'on';
end

% Initial settings
updateAxesWithEngUnits(hfig);
updateStatus(hfig);
centerCursors([],[],hfig);

end

% ---------------------------------------------------------------
% H E L P    S Y S T E M
% --------------------------------------------------------------

%--------------------------------------------------------------
function installContextHelp(hfig)

ud = get(hfig,'UserData');

main = {'Label','&What''s This?', 'Callback',@HelpGeneral, 'Parent'};

setWTC(hfig, main, [ud.hsurface ud.hax(1)], 'Spectrogram image');
setWTC(hfig, main, ud.hthumb, 'Zoom Window Panner');
setWTC(hfig, main, [ud.himageCbar ud.hax(3)], 'Colorbar');
setWTC(hfig, main, ud.htextStatus, 'Status Bar');
setWTC(hfig, main, ud.htextAmp, 'Magnitude Readout');
setWTC(hfig, main, ud.htextFreq, 'Frequency Readout');
setWTC(hfig, main, ud.htextTime, 'Time Readout');
setWTC(hfig, main, [ud.htimePlot ud.hax(2)], 'Time Plot');

setWTC(hfig,main, ud.hspecX, 'Frequency Crosshair');
setWTC(hfig,main, [ud.htimeY ud.hspecY], 'Time Crosshair');
setWTC(hfig,main, ud.hfig, 'Spectrogram Demo');

end

%--------------------------------------------------------------
function setWTC(hfig,main,hItem,tagStr) %#ok<INUSD,INUSL>
% setWT Set the "What's This?" context menu and callback:
hc = uicontextmenu('Parent',hfig);
% Do not add context help for now, we will support this later.
%uimenu(main{:},hc, 'Tag',['WT?' tagStr]);
set(hItem,'UIContextMenu',hc);

end

% ---------------------------------------------------------------
% C O N T E X T   M E N U S
% --------------------------------------------------------------

%-----------------------------------------------------------------
function installContextMenus(hfig)

installSpecgramModeMenus(hfig);
installColorbarMenus(hfig);
installTimePannerMenus(hfig);

end

%-----------------------------------------------------------------
function installSpecgramModeMenus(hfig)

% Additional menus to prepend to the spectrogram context menu:

ud = get(hfig,'UserData');
hc = get(ud.hsurface,'UIContextMenu');

% Update the menu on-the-fly:
set(hc,'Callback', @specgramMenuRenderCallback);

hEntry=[];

opts={hc,getStrFromCat('fullView'),@zoomFull, 'Tag', 'fullViewCtxMenu'};
hEntry(end+1) = createContext(opts);
opts={hc,getStrFromCat('centerCursors'),@centerCursors, 'Separator','on',...
  'Tag', 'centerCursorstxMenu'};
hEntry(end+1) = createContext(opts);
opts={hc,getStrFromCat('cascadePlot'),@toggle2D3DView, 'Separator','on',...
  'Checked','off','Tag', 'cascadePlotCtxMenu'};
hEntry(end+1) = createContext(opts);
opts={hc,getStrFromCat('rotateLeft'),@rotate3DViewLeft,...
  'Tag', 'rotLeftCtxMenu'};
hEntry(end+1) = createContext(opts);
opts={hc,getStrFromCat('rotateRight'),@rotate3DViewRight,...
  'Tag', 'rotRightCtxMenu'};
hEntry(end+1) = createContext(opts);
opts={hc,getStrFromCat('reset3DView'),@reset3DView,...
  'Tag', 'reset3DViewCtxMenu'};
hEntry(end+1) = createContext(opts);

% Give each menu item a vector of handles to all peer menus
set(hEntry,'UserData',hEntry);

end

%-----------------------------------------------------------------
function installColorbarMenus(hfig)
% Additional menus to prepend to the colorbar context menu:

ud = get(hfig,'UserData');
hc = get(ud.himageCbar,'UIContextMenu');  % ud.hax(1) also?

opts={hc,getStrFromCat('colormap'),''};
hCmap = createContext(opts);

hEntry=[];  % holds handles to each colormap menu item
opts={hCmap,'Parula',@changeCMap, 'Checked','on'};
hEntry(end+1) = createContext(opts);
opts={hCmap,'Jet',@changeCMap};
hEntry(end+1) = createContext(opts);
opts={hCmap,'Hot',@changeCMap};
hEntry(end+1) = createContext(opts);
opts={hCmap,'Gray',@changeCMap};
hEntry(end+1) = createContext(opts);
opts={hCmap,'Bone',@changeCMap};
hEntry(end+1) = createContext(opts);
opts={hCmap,'Copper',@changeCMap};
hEntry(end+1) = createContext(opts);
opts={hCmap,'Pink',@changeCMap};
hEntry(end+1) = createContext(opts);

opts={hc,getStrFromCat('setLimits'),@manualCmapLimits, 'Separator','on'};
createContext(opts);

opts={hc,getStrFromCat('resetLimits'),@resetCmapLimits};
createContext(opts);

% Give each menu item a vector of handles to all peer menus
set(hEntry,'UserData',hEntry);

%fixupContextOrder(hc);

end

%-----------------------------------------------------------------
function installTimePannerMenus(hfig)

% Additional menus to prepend to the time-panner context menu:

ud = get(hfig,'UserData');
hthumb = ud.hthumb;  % add to time axis as well?
hc = get(hthumb, 'UIContextMenu');

% Update the menu on-the-fly:
set(hc,'Callback', @focusMenuRenderCallback);

hEntry=[];  % holds handles to each colormap menu item

opts={hc,getStrFromCat('focusIn'),@focusTimeIn};
hEntry(end+1) = createContext(opts);

opts={hc,getStrFromCat('previousFocus'),@focusTimePrev};
hEntry(end+1) = createContext(opts);

opts={hc,getStrFromCat('resetFocus'),@focusTimeReset};
hEntry(end+1) = createContext(opts);

% Give each menu item a vector of handles to all peer menus
set(hEntry,'UserData',hEntry);

%fixupContextOrder(hc);

updateFocusHistoryMenu(hfig); % pass any focus context menu

end

%-----------------------------------------------------------------
function specgramMenuRenderCallback(~, ~)

hfig=gcbf; hparent=gcbo;
ud = get(hfig,'UserData');

hAllMenus = get(hparent,'Children'); % vector of handles to context menus

tags = {hAllMenus.Tag};
hRotateRight = hAllMenus(strcmp(tags, 'rotRightCtxMenu'));
hRotateLeft = hAllMenus(strcmp(tags, 'rotLeftCtxMenu'));
hReset = hAllMenus(strcmp(tags, 'reset3DViewCtxMenu'));

if ud.View3DStatus
  ena = 'on';
else
  ena = 'off';
end
set(hRotateRight,'Enable',ena);
set(hRotateLeft,'Enable',ena);
set(hReset,'Enable',ena);
end

%-----------------------------------------------------------------
function hMenu = createContext(opts)
% Helper function to append additional context menus
args = {'Parent',opts{1}, 'Tag',opts{2}, 'Label',opts{2},'Callback',opts{3:end}};
hMenu = uimenu(args{:});
end

%-----------------------------------------------------------------
function fixupContextOrder(hContext) %#ok<DEFNU>
% Put the first context menu entry (the "What's This?" entry)
%  last in the context menu list, and turn on the separator
%  for the "What's This?" entry
childList = get(hContext,'Children');
childList = childList([end 1:end-1]);
set(hContext,'Children',childList);
set(childList(1),'Separator','on');

end

%---------------------------------------------------------------
function changeCMap(~,~)

hco=gcbo; hfig=gcbf;
% Reset checks on all colormap menu items:
set(get(hco,'UserData'),'Checked','off');
set(hco,'Checked','on');

% Update figure colormap:
cmapStr = lower(get(hco,'Label'));
cmap = feval(cmapStr);
set(hfig,'Colormap',cmap);

end

% ---------------------------------------------------------------
% F O C U S    S Y S T E M
% --------------------------------------------------------------

%---------------------------------------------------------------
function pushCurrToFocusHistory(hfig)

ud = get(hfig,'UserData');
haxTime = ud.hax(2);

% focus history is stored in userdata of time-panner axis as either an
% empty vector or cell, or as a cell-array of 2-element x-lim vector.

% get current time-axis limits
currXlim = get(haxTime,'XLim');

currHistory = get(haxTime,'UserData');
if isempty(currHistory)
  updated_focus_history = {currXlim};
else
  updated_focus_history = [currHistory {currXlim}];
end
set(haxTime,'UserData',updated_focus_history);

updateFocusHistoryMenu(hfig);

end

%---------------------------------------------------------------
function hist_xlim = popFromFocusHistory(hfig)

ud = get(hfig,'UserData');
haxTime = ud.hax(2);
currXlim = get(haxTime,'XLim'); % get current time-axis limits

currHistory = get(haxTime,'UserData');
if isempty(currHistory)
  % no prev focus info recorded
  warning(message('signal:rpmmapplot:popEmptyFocusStack'));
  hist_xlim = currXlim;
  
else
  % Pop last history xlim
  hist_xlim = currHistory{end};
  currHistory(end) = [];
  set(haxTime,'UserData',currHistory);
end

updateFocusHistoryMenu(hfig);

end

%---------------------------------------------------------------
function clearFocusHistory(hfig)
% Remove all previous focus entries

ud = get(hfig,'UserData');
haxTime = ud.hax(2);
set(haxTime,'UserData',[]);

updateFocusHistoryMenu(hfig);

end

%---------------------------------------------------------------
function updateFocusHistoryMenu(hfig)

ud = get(hfig,'UserData');
haxTime = ud.hax(2);

% Update 'Previous Focus' context menu label
currHistory = get(haxTime,'UserData');
histLen = length(currHistory);
str = getStrFromCat('prevFocus');
if histLen>0
  str = [str ' (' num2str(histLen) ')'];
  ena = 'on';
else
  ena = 'off';
end

% Get panner context menu handle:
hmenu = findobj( get(get(ud.hthumb, 'UIContextMenu'),'Children'),'Tag',getStrFromCat('focusIn'));
hAllMenus = get(hmenu,'UserData'); % vector of handles to context menus
hFocusPrev = hAllMenus(2);
set(hFocusPrev, 'Label',str);
set(hAllMenus(2:3), 'Enable',ena);  % Prev and Reset Focus menus

end 

%---------------------------------------------------------------
function focusMenuRenderCallback(~, ~)
% Used to update the enable of the "Focus In" menu item. Only enabled if
% thumb_xlim ~= curr_xlim

hfig=gcbf; hparent=gcbo;
ud = get(hfig,'UserData');
hAllMenus = get(hparent,'Children'); % vector of handles to context menus

% Enable 'Focus on Window' if zoom window is less than entire panner
%
hFocusIn = hAllMenus(end);  % 'Focus on Zoom' entry
haxTime = ud.hax(2);
currXlim = get(haxTime,'XLim'); % get current time-axis limits
% Get thumbnail xlim vector:
thumbXdata = get(ud.hthumb,'XData');  % current thumbnail patch coords
thumbXlim  = [min(thumbXdata) max(thumbXdata)]; % convert to xlim
if ~isequal(currXlim, thumbXlim)
  ena='on';
else
  ena='off';
end
set(hFocusIn,'Enable',ena);

end 

%---------------------------------------------------------------
function focusTimeIn(~,~)

hfig=gcbf;

% get current time-axis (panner) limits
ud = get(hfig,'UserData');
haxTime = ud.hax(2);
currXlim = get(haxTime,'XLim');

% Get thumbnail xlim vector:
thumbXdata = get(ud.hthumb,'XData');  % current thumbnail patch coords
thumbXlim  = [min(thumbXdata) max(thumbXdata)]; % convert to xlim

if ~isequal(currXlim, thumbXlim)
  pushCurrToFocusHistory(hfig);
  
  % Zoom in to thumb limits
  haxTime = ud.hax(2);
  
  set(haxTime,'XLim', thumbXlim);
  updateAxesWithEngUnits(gcbf);
end

end

%---------------------------------------------------------------
function focusTimePrev(~,~)

hfig=gcbf;
ud = get(hfig,'UserData');
haxTime = ud.hax(2);

% Reset to last focus
xlim = popFromFocusHistory(hfig);

set(haxTime, 'XLim',xlim);
updateAxesWithEngUnits(gcbf);

end 

%---------------------------------------------------------------
function focusTimeReset(~,~,hfig)
% Remove all previous focus entries

if nargin < 3
  hfig=gcbf;
end
clearFocusHistory(hfig);

% Reset focus zoom:
ud = get(hfig,'UserData');
hax_time = ud.hax(2);
surfXdata = get(ud.hsurface,'XData');

if length(surfXdata) == 1
  set(hax_time, 'XLim', [0 surfXdata]);
else
  set(hax_time,'XLim',[min(surfXdata), max(surfXdata)]);
end
updateAxesWithEngUnits(hfig);

end 

% ---------------------------------------------------------------
% UTILITY FUNCTIONS
% --------------------------------------------------------------
function  hFrame = createUIFramework(res, type)
%createUIFramework Create the uiframework.

% Create uimgr UI frame
hFrame = CreateUIFrame(type);

CreateBaseMenus(hFrame);
CreateBaseToolbar(hFrame);
ud = CreateBaseStatusbar(hFrame);
% Render ui objets

% Top of the uimgr tree
% technically, this is the only handle we really need to keep
% all others could be found from this, using hFrame.findchild(...)
ud.huiframe = hFrame;
% Other params in the userdata:
ud.bGraphic = false;
ud.param_dlg = [];
if strcmpi(type,'Order')
  ud.ResUnits = getStrFromCat('resReadoutUnits');
else
  ud.ResUnits = 'Hz';
end
ud.ResTxt = [getStrFromCat('resReadoutLbl') ': '];
ud.ResValue = sprintf('%.3f',res);
ud.MapType = type;

ud.View3DStatus = false;

set(hFrame,...
  'Tag', 'rpmmapFigTag',...
  'CloseRequestFcn', @closeCb, ...
  'UserData', ud);

hPrintBehavior = hggetbehavior(hFrame,'Print');
set(hPrintBehavior,'WarnOnCustomResizeFcn','off');

end
%----------------------------------------------------------
function CreateBaseMenus(hFrame)

% Menus group
%hm = uimgr.uimenugroup('Menus');
%fig = figure;
%hm = uimenu(hFrame,'Label','Menus');

% Files
mFile = uimenu(hFrame,'Label',getStrFromCat('menuFile'),...
        'Tag','uimgr.uimenugroup_File');
uimenu(mFile,'Label', getStrFromCat('menuPrint'),...
      'Tag','uimgr.uimenugroup_Print','Callback',@printdlgCb);
uimenu(mFile,'Label', getStrFromCat('menuPrintPreview'),...
       'Tag','uimgr.uimenugroup_Printview','Callback',@printpreviewCb);
                
uimenu(mFile,'Text', getStrFromCat('menuClose'),...
        'Separator','on','Tag','uimgr.uimenugroup_FileOpt',...
        'Callback', @closeCb);

% Tools
mTool = uimenu(hFrame,'Label',getStrFromCat('menuTools'));

uimenu(mTool,'Label', getStrFromCat('menuChangeTitle'),...
                'Tag','changeTitleToolsMenuItem',...
                'Callback', @changeTitle);
%mTitlegroup = uimenugroup(mTool,'Titlegroup',mTitle);

uimenu(mTool,'Label', getStrFromCat('menuFullView'),...
                   'Tag','fullViewToolsMenuItem',...
                   'Separator','on',...
                   'Callback', @zoomFull);
uimenu(mTool,'Label', getStrFromCat('menuZoomIn'),...
       'Tag','zoomInToolsMenuItem',...
       'Callback', @toggleZoom);
uimenu(mTool,'Label', getStrFromCat('menuZoomOut'),...
       'Tag','zoomOutToolsMenuItem',...
       'Callback', @toggleZoom);

uimenu(mTool,'Label', getStrFromCat('menuZoomX'),...
    'Tag','zoomXToolsMenuItem','Callback', @toggleZoom);

uimenu(mTool,'Label', getStrFromCat('menuZoomY'),...
       'Callback', @toggleZoom,'Tag','zoomYToolsMenuItem');

uimenu(mTool,'Label', getStrFromCat('menuCenterCursors'),...
       'Tag','centerCursorsToolsMenuItem',...
       'Separator','on',...
       'Callback', @centerCursors);
uimenu(mTool, 'Label', getStrFromCat('menuCascadePlot'),...
       'Checked', 'off', 'Tag','cascadePlotToolsMenuItem',...
       'Separator','on',...
       'Callback', @toggle2D3DView);
uimenu(mTool, 'Label', getStrFromCat('rotateLeft'),...
       'Tag','rotLeftToolsMenuItem','Enable','off',...
        'Callback', @rotate3DViewLeft);
uimenu(mTool, 'Label', getStrFromCat('rotateRight'),...
       'Tag','rotRightToolsMenuItem','Enable','off',...
       'Callback', @rotate3DViewRight);
uimenu(mTool, 'Label', getStrFromCat('reset3DView'),...
       'Callback', @reset3DView,...
       'Tag','reset3DViewToolsMenuItem',...
       'Enable','off');
end

%----------------------------------------------------------
function CreateBaseToolbar(hFrame)

icon_file = 'rpmmapploticons.mat';
icons = spcwidgets.LoadIconFiles(icon_file);

ht = uitoolbar(hFrame);

uipushtool(ht,'CData', icons.fit_to_view,...
              'TooltipString',getStrFromCat('toolbarFullView'),...
              'ClickedCallback', @zoomFull,...
              'Tag', 'fullViewToolbarButton');
uitoggletool(ht,'CData', icons.zoom_in,...
            'Tag','zoomInToolbarButton','TooltipString',getStrFromCat('toolbarZoomIn'), ...
            'ClickedCallback',@toggleZoom );

uitoggletool(ht,'CData', icons.zoom_out,...
            'Tag','zoomOutToolbarButton','TooltipString',getStrFromCat('toolbarZoomOut'), ...
            'ClickedCallback',@toggleZoom );

uitoggletool(ht,'CData', icons.toggle_zoom_x,...
            'Tag','zoomXToolbarButton','TooltipString',getStrFromCat('toolbarZoomX'), ...
            'ClickedCallback',@toggleZoom );

uitoggletool(ht,'CData', icons.toggle_zoom_y,...
            'Tag','zoomYToolbarButton','TooltipString',getStrFromCat('toolbarZoomY'), ...            
            'ClickedCallback',@toggleZoom );

uipushtool(ht,'Tag','centerCursorsToolbarButton',...
           'CData',icons.center_crosshair,...
           'TooltipString',getStrFromCat('toolbarCenterCursors'),...      
           'Separator', 'on', ...
           'ClickedCallback',@centerCursors);
        
uitoggletool(ht,'CData', icons.waterfall,...
            'Tag','cascadePlotToolbarButton','TooltipString',getStrFromCat('toolbarCascadePlot'), ...
            'Separator', 'on', ...
            'ClickedCallback',@toggle2D3DView );

        
uipushtool(ht,'CData', icons.rotateLeft,...
            'Tag','rotLeftToolbarButton','TooltipString',getStrFromCat('rotateLeft'), ...
            'ClickedCallback',@rotate3DViewLeft,'Enable','off');

uipushtool(ht,'CData', icons.rotateRight,...
            'Tag','rotRightToolbarButton','TooltipString',getStrFromCat('rotateRight'), ...
            'ClickedCallback',@rotate3DViewRight,'Enable','off');

uipushtool(ht,'CData', icons.rotateReset,...
            'Tag','resetViewToolbarButton','TooltipString',getStrFromCat('reset3DView'), ...
            'ClickedCallback',@reset3DView,'Enable','off');

end

%----------------------------------------------------------
function ud = CreateBaseStatusbar(hFrame)

hMainFlow = uiflowcontainer('v0', ...
    'Parent', hFrame, ...
    'Tag', 'ApplicationLayoutContainer', ...
    'Flowdirection', 'topdown', ...
    'HitTest', 'off', ...
    'Margin', 0.1);

% Insert container for visualization area
uicontainer('Parent', hMainFlow, ...
    'Tag', 'VisualizationAreaContainer', ...
    'HitTest', 'off');

grayLine = uicontrol('Parent', hMainFlow, ...
    'Style', 'frame', ...
    'ForegroundColor', get(0,'FactoryUipanelShadowColor'), ...
    'Tag', 'GraySeparatorLine', ...
    'Position', [0 0 1 1]);
set(grayLine, 'HeightLimits', [0 0]);
whiteLine = uicontrol('Parent', hMainFlow, ...
    'Style', 'frame', ...
    'ForegroundColor', get(0,'FactoryUipanelHighlightColor'), ...
    'Tag', 'WhiteSeparatorLine', ...
    'Position', [0 0 1 1]);
set(whiteLine, 'HeightLimits', [0 0]);

setappdata(hMainFlow, 'LineSeparators', [whiteLine, grayLine]);
%hs = uimgr.uistatusbar('StatusBar');
            
hs = spcwidgets.StatusBar(hMainFlow,'Text','',...
     'Tag','StatusBar');
ud.hStatusBar = hs;
ud.hStatus = spcwidgets.Status(hs,'Text','','Tag','Res',...
           'Tooltip', getStrFromCat('resReadoutToolTip'), ...
           'Width', 160);
% Use these if we want cursor readouts in status bar:
% hsTime = uimgr.uistatus('Time');
% hsTime.setWidgetPropertyDefault(...
%   'Text', 'T: ', ...
%   'Tooltip', 'Time at vertical cursor.', ...
%   'Width', 100);
% hs.add(hsTime);
%
% hsFreq = uimgr.uistatus('Freq');
% hsFreq.setWidgetPropertyDefault(...
%   'Text', 'F: ', ...
%   'Tooltip', 'Frequency at horizontal cursor.', ...
%   'Width', 100);
% hs.add(hsFreq);
%
% hsRPM = uimgr.uistatus('RPM');
% hsRPM.setWidgetPropertyDefault(...
%   'Text', 'RPM: ', ...
%   'Tooltip', 'RPM at vertical cursor.', ...
%   'Width', 100);
% hs.add(hsRPM);
%
% hsAmp = uimgr.uistatus('Amp');
% hsAmp.setWidgetPropertyDefault(...
%   'Text', 'A: ', ...
%   'Tooltip', 'Amplitude at cursors intersection.', ...
%   'Width', 100);
% hs.add(hsAmp);

end

%----------------------------------------------------------
function hFrame = CreateUIFrame(type)


if strcmpi(type,'Order')
    titStr = getStrFromCat('orderMap');
else
    titStr = getStrFromCat('freqMap');
end

hFrame = figure('Tag', 'uimgr.uifigure_UIFrame',...
 'MenuBar','none','ToolBar','none',...
 'NumberTitle', 'off','Name',titStr,...
 'ResizeFcn',@resizeFig, ...
 'Position',[624, 236, 720, 540],...
 'PaperPositionMode','auto',...
 'DockControls','Off','Visible','off');

end

%----------------------------------------------------------
function update3DControlItems(hfig)

ud = get(hfig,'UserData');
flag = ud.View3DStatus;

tags = {...  
  'rotLeftToolsMenuItem',...
  'rotRightToolsMenuItem',...
  'reset3DViewToolsMenuItem',...
  'rotLeftCtxMenu',...  
  'rotRightCtxMenu',...
  'reset3DViewCtxMenu',...
  'rotLeftToolbarButton',...
  'rotRightToolbarButton',...
  'resetViewToolbarButton'};
  
if flag
  % We are in 3D mode
  % Enable tool menu items and check the cascade plot checkbox
  enab = 'on';  
else
  enab = 'off'; 
end

item = findobj('Tag','cascadePlotToolsMenuItem');
set(item,'Checked',enab);
item = findobj('Tag','cascadePlotCtxMenu');
set(item,'Checked',enab);
item = findobj('Tag','cascadePlotToolbarButton');
set(item,'State',enab);

for idx = 1:numel(tags)
  item = findobj('Tag',tags{idx});
  set(item,'Enable',enab);
end

end
%----------------------------------------------------------
function updateZoomControlItems(hfig)

% Check if zoom is on, and the zoom state and synchronize all menu and
% toolbar controls
ud = get(hfig,'UserData');
zoomOn = ud.ZoomOn;
zoomMode = ud.ZoomMode;

tagsZoomIn  = {'zoomInToolsMenuItem','zoomInToolbarButton','zoomInCtxMenu'};
tagsZoomOut = {'zoomOutToolsMenuItem','zoomOutToolbarButton','zoomOutCtxMenu'};
tagsZoomX   = {'zoomXToolsMenuItem','zoomXToolbarButton','zoomXCtxMenu'};
tagsZoomY   = {'zoomYToolsMenuItem','zoomYToolbarButton','zoomYCtxMenu'};

if zoomOn
  switch zoomMode    
    case 'ZoomIn'
      setCtrlState(hfig,tagsZoomIn,'on');
      setCtrlState(hfig,[tagsZoomOut,tagsZoomX,tagsZoomY],'off');
    case 'ZoomOut'
      setCtrlState(hfig,tagsZoomOut,'on');
      setCtrlState(hfig,[tagsZoomIn,tagsZoomX,tagsZoomY],'off');
    case 'ZoomX'
      setCtrlState(hfig,tagsZoomX,'on');
      setCtrlState(hfig,[tagsZoomIn,tagsZoomOut,tagsZoomY],'off');
    case 'ZoomY'
      setCtrlState(hfig,tagsZoomY,'on');
      setCtrlState(hfig,[tagsZoomIn,tagsZoomOut,tagsZoomX],'off'); 
  end
else
  setCtrlState(hfig,[tagsZoomIn,tagsZoomOut,tagsZoomX,tagsZoomY],'off');
end
end
%----------------------------------------------------------
function status = getZoomModeFromTag(ctrl)

switch ctrl.Tag
  case {'zoomInToolsMenuItem','zoomInToolbarButton','zoomInCtxMenu'}
    status = 'ZoomIn';          
  case {'zoomOutToolsMenuItem','zoomOutToolbarButton','zoomOutCtxMenu'}
    status = 'ZoomOut';
  case {'zoomXToolsMenuItem','zoomXToolbarButton','zoomXCtxMenu'}
    status = 'ZoomX';
  case {'zoomYToolsMenuItem','zoomYToolbarButton','zoomYCtxMenu'}
    status = 'ZoomY';
  otherwise
    status = [];
end
end

%----------------------------------------------------------
function setCtrlState(hfig,tags,state)
% Set checked or state property of control with tag in tags cell array to
% the state in state input.

for idx = 1:numel(tags)
  
  item = findobj(hfig,'Tag',tags{idx});
  if ~isempty(item)
        
    if isprop(item,'Checked')
      item.Checked = state;
    elseif isprop(item,'State')
      item.State = state;
    end
    
  end
  
end
end
%----------------------------------------------------------
function changeTitle(~,~,hfig)
if nargin < 3
  hfig = gcbf;
end

ud = get(hfig,'UserData');
haxSpec = ud.hax(1);

dlgTitle = getStrFromCat('title');
prompt={[getStrFromCat('title') ':']};
tit = get(haxSpec,'Title');
def = {tit.String};
numLines = 1;
strs = inputdlg(prompt,dlgTitle,numLines,def);

if isempty(strs)
  return
end
tit.String = strs{1};

end
%----------------------------------------------------------
  function str = getStrFromCat(inStr)
    % Get string from catalog
        
    str = getString(message(['signal:rpmmapplot:' inStr]));
    
  end
  
  %----------------------------------------------------------
  function m = getMultiplier(units)
  
  if isempty(units)
    m = '';
    return;
  end
  
  switch units
    case 'y'
      m = '\times 1e-24';
    case 'z'
      m = '\times 1e-21';
    case 'a'
      m = '\times 1e-18';
    case 'f'
      m = '\times 1e-15';    
    case 'p'
      m = '\times 1e-12';
    case 'n'
      m = '\times 1e-9';
    case '\mu'
      m = '\times 1e-6';
    case 'm'
      m = '\times 1e-3';

    case 'k'
      m = '\times 1e3';
    case 'M'
      m = '\times 1e6';
    case 'G'
      m = '\times 1e9';
    case 'T'
      m = '\times 1e12';
    case 'P'
      m = '\times 1e15';
    case 'E'
      m = '\times 1e18';
    case 'Z'
      m = '\times 1e21';
    case 'Y'
      m = '\times 1e24';
  end
  end
  

