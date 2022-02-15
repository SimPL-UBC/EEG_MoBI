function [varargout] = rpmmap(xIn,fsIn,rpmIn,maptype,varargin)
%RPMMAP Compute frequency or order maps vs. rpm values.

% Copyright 2015-2019 The MathWorks, Inc.
%#codegen

% parse input parameters
inpArgs = cell(size(varargin));
if nargin > 4
    [inpArgs{:}] = convertStringsToChars(varargin{:});
else
    inpArgs = varargin;
end

[resol,winName,winparam,amplitude,scale,overlapPercent] = ...
  parseOptions(maptype,xIn,fsIn,rpmIn,inpArgs{:});

% cast to enforce precision rules (we already checked that the inputs are
% numeric).
rpm = double(reshape(rpmIn,[],1));
fs  = double(fsIn(1));
x   = reshape(xIn,[],1);
% Define time vector
time = ((0:length(x)-1)/fs).';
if strcmp(maptype,'order')
  % Convert signal to rotational or order domain.
  % xp is the resampled signal (constant samples/cycle) and fsp is the
  % samples/cycle rate. phaseUp, rpmUp and timeUp vectors are at an
  % upsampled rate of 15*fs.
  [xp, fsp, phaseUp, rpmUp, timeUp] = toConstantSamplesPerCycle(x,fs,rpm,time);
else
  xp = x;
  fsp = fs;
end

% Compute DFT window length
% Define the minimum and maximum allowed window lengths
minWindowLength = 4;
maxWindowLength = length(xp);

if isempty(resol)
  % Use the default resolution value: (sampling frequency)/128 for
  % frequency maps and (sampling frequency)/256 for order maps.
  if strcmp(maptype,'order')
    resolution = fsp/256;
  else
    resolution = fsp/128;
  end
  % Make sure resolution is inside range of allowed values
  % If it is not, use the minimum or maximum allowed window length
  [isResolValid,winLengthBound] = ...
  validateResolutionValue(resolution,minWindowLength,maxWindowLength,fsp,winName,winparam,false);
else
   resolution = resol;
  % Validate input resolution value
  [isResolValid,winLengthBound] = ...
  validateResolutionValue(resolution,minWindowLength,maxWindowLength,fsp,winName,winparam,true);
end

if isResolValid
  % We have a valid resolution so compute window length
  [~,winLength] = signal.internal.getWinDurationForAGivenRBW(...
    resolution,winName,winparam,fsp,true);
else
  % We have an invalid resolution so use the upper or lower bound
  winLength = winLengthBound;
end

% Compute the number of overlapping samples 
nOverlap = min(ceil(overlapPercent/100*winLength),winLength-1);

% Create DFT window and normalize it
win = getWindow(winName, winLength,winparam);
win = win/sum(win);

% Compute the resolution of the window and the FFT length
winres = enbw(win)*fsp/winLength;
nfft = max(256,length(win));

% Generate the frequency/order map
if strcmp(maptype,'order')
  [mapUnscaled,ordSpec,phaseSpec] = spectrogram(xp,win,nOverlap,nfft,fsp,'onesided');
  phaseSpec = phaseSpec(:);
else
  [mapUnscaled,freqSpec,timeSpec] = spectrogram(xp,win,nOverlap,nfft,fsp,'onesided');
  timeSpec = timeSpec(:);
end

% Apply amplitude and scale of the map
map = mapAmplitudeScale(mapUnscaled,amplitude,scale,nfft);
    
if strcmp(maptype,'order')
  % Limit to max order value as we oversampled in the order domain for
  % better results. Recall that Omax = fs/(2*max(rpm/60))
  idxTemp = find(ordSpec <= fs/(2*max(rpm/60)));
  if idxTemp(end)+1 <= length(ordSpec)
    % Get one extra order point to ensure we have the max order value
    idx = [idxTemp;idxTemp(end)+1];
  else
    idx = idxTemp;
  end
  ordSpec = ordSpec(idx);
  map = map(idx,:);
end

%---------------------------------------------------------------------
% Plot or assign outputs
%---------------------------------------------------------------------
% Create plot of ordermap if no output arguments are specified
% Otherwise, assign output variables to suppress command window output
if nargout == 0
  coder.internal.assert(coder.target('MATLAB'),'signal:rpmmap:PlottingNotSupported');
  % Surf encodes data points at the edges and takes the color from the last
  % edge so we need to add an additional point so that surf does the right
  % thing. This is important especially when spectrogram has only one
  % estimate (e.g. window length = signal length). For the plot we set time
  % or phase values to be at: nwin/2-a/2, nwin/2+a/2, nwin/2+3a/2,
  % nwin/2+5a/2 ... where a is the number of new samples for each segment
  % (i.e. nwin-noverlap). For the case of zero overlap this corresponds to
  % 0, nwin, 2*nwin, ...
  a = winLength - nOverlap;
  if strcmpi(maptype,'order')
    phaseVectForPlot = [(winLength/2-a/2)/fsp;  phaseSpec + ((a/2)/fsp)];
    rpmVectForPlot   = interp1(phaseUp, rpmUp, phaseVectForPlot,'linear','extrap');
    timeVectForPlot  = interp1(phaseUp, timeUp, phaseVectForPlot,'linear','extrap');
    map              = [map map(:,end)];
    
    rpmmapplot(map, ordSpec, timeVectForPlot, rpmVectForPlot, ...
      winres,'order', scale, amplitude, rpm, time);
  else
    timeVectForPlot = [(winLength/2-a/2)/fsp; timeSpec + ((a/2)/fsp)];
    rpmVectForPlot  = interp1(time, rpm, timeVectForPlot,'linear','extrap');
    map             = [map map(:,end)];
    
    rpmmapplot(map, freqSpec, timeVectForPlot, rpmVectForPlot, ...
      winres,'frequency', scale, amplitude, rpm, time);
    
  end
  
else
  
  % Assign output variables
  varargout{1} = map;
  
  if strcmp(maptype,'order')
    % Order map outputs are [map, order, rpm, time, res]
    if nargout > 1
      varargout{2} = ordSpec;
      if nargout > 2
        % Interpolate RPM values based on upsampled phase and rpm values,
        % and on spectrogram phase output vector
        rpmSpecDerived = interp1(phaseUp,rpmUp,phaseSpec,'linear','extrap');
        varargout{3} = rpmSpecDerived;
        if nargout > 3
          % Interpolate time values based on upsampled phase and rpm values,
          % and on spectrogram phase output vector
          timeSpecDerived = interp1(phaseUp,timeUp,phaseSpec,'linear','extrap');
          varargout{4} = timeSpecDerived;
          if nargout > 4
            varargout{5} = winres;
          end
        end
      end
    end
  else
    if nargout > 1
      varargout{2} = freqSpec;
      if nargout > 2
        % Interpolate RPM values based on time and rpm values, and on
        % spectrogram time output vector
        rpmSpecDerived = interp1(time,rpm,timeSpec,'linear','extrap');
        varargout{3} = rpmSpecDerived;
        if nargout > 3
          varargout{4} = timeSpec;
          if nargout > 4
            varargout{5} = winres;
          end
        end
      end
    end
  end
end
end

%--------------------------------------------------------------------------
function [xp, fsp, phaseUp, rpmUp, timeUp] = toConstantSamplesPerCycle(x,fs,rpm,time)

% Compute the maximum order that can be present in X with no aliasing. Max
% frequency, fmax, is the signal's max frequency (max(rpm/60)) multiplied
% by max order. For Nyquist to hold fs > 2*fmax = 2*max(rpm/60)*Omax.
Omax = fs/(2*max(rpm/60));

% Define sampling rate (samples/cycle) in phase domain based on the maximum
% order. Make sampling rate 4 times the Nyquist rate in the order domain.
fsp = 4*(2*Omax);

% Define upsample factor. Upsampling will improve accuracy when converting
% from constant samples/second to constant samples/cycle. In the worst case
% scenario, when time signal is critically sampled in time at Fmax/2, we
% are increasing the Nyquist frequency by 15.
upFactor = 15;

% Upsample x and rpm
if isa(x,'single')
  xtemp = resample(double(x),upFactor,1);
  xUp = single(xtemp);
else
  xUp = resample(x,upFactor,1);
end
% Get upsampled time and rpm vectors
timeUp = (0:length(xUp)-1).'/(upFactor*fs);
rpmUp = interp1(time, rpm, timeUp, 'linear','extrap');

% Estimate the phase of each signal sample by integrating the instantaneous
% signal frequency which is rpmUp/60. Divide by sampling rate which is
% upFactor*fs;
phaseUp = cumtrapz(rpmUp/(60*upFactor*fs));

% Interpolate signal x at constant phase increments (i.e. constant
% samples/cycle). xp is uniformly sampled in the rotational domain --> same
% samples per rotation for any rpm. Use only strictly increasing values of
% phaseUp (values may not be unique for very small or very large rpm
% inputs).
constPhase = (phaseUp(1):1/fsp:(phaseUp(end)))';
strictincrIdx = [true; diff(phaseUp) > 0];
phaseUp = phaseUp(strictincrIdx);
rpmUp = rpmUp(strictincrIdx);
timeUp = timeUp(strictincrIdx);
xp = interp1(phaseUp, xUp(strictincrIdx), constPhase, 'linear','extrap');
end

%--------------------------------------------------------------------------
function [resolution, win, winparam, amplitude, scale, overlapPercent, funName] = parseOptions(maptype,x,fs,rpm,varargin)

% default values for n-v pairs
if strcmp(maptype,'order')
    defaultWin = 'flattopwin';
    funName = 'rpmordermap';
else
    defaultWin = 'hann';
    funName = 'rpmfreqmap';
end
defaultAmp = 'rms';
defaultScale= 'linear';
defaultOverlapPercent = 50;
% Check valid values for x,fs,rpm
validateattributes(x,{'single','double'},...
    {'real','nonsparse','finite','vector'},funName,'X');
validateattributes(fs,{'numeric'},...
    {'real','positive','nonsparse','finite','scalar'},funName,'Fs');
validateattributes(rpm,{'numeric'},...
    {'real','positive', 'nonsparse','finite','vector'},funName,'RPM');

% Check valid input dimensions
coder.internal.assert(length(x) >= 18,'signal:rpmmap:MustBeMinLength');

%check rpm is the same size and x
coder.internal.assert(length(x) == length(rpm),'signal:rpmmap:MustBeSameLength');

%check if a resolution parameter is specified and check if valid value
if ~isempty(varargin) && isnumeric(varargin{1})
    validateattributes(varargin{1},{'numeric'},...
        {'real','positive','nonsparse','finite','scalar'},funName,'resolution');
    resolution = double(varargin{1}(1));
    args = {varargin{2:end}};
else
    args = varargin;
    resolution = [];
end

if ~isempty(args)
    % args should have even number of elements. 
    coder.internal.assert(~isodd(numel(args)),'signal:rpmmap:NVMustBeEven');
    % Parse the argument list to obtain the name-value pairs
    if coder.target('MATLAB')
        p = inputParser;
        p.FunctionName = funName;
        addParameter(p,'Amplitude',defaultAmp);
        addParameter(p,'Scale',defaultScale);
        addParameter(p,'OverlapPercent',defaultOverlapPercent);
        addParameter(p,'Window',defaultWin);
        parse(p,args{:});
        ampOut     =  p.Results.Amplitude;
        scaleOut   =  p.Results.Scale;
        ovpOut     =  p.Results.OverlapPercent;
        winOut     =  p.Results.Window; % winOut can be a string or a cell Array       
    else
        params = struct('Amplitude',uint32(0),'Scale',uint32(0),...
                         'OverlapPercent',uint32(0),'Window',uint32(0));
        poptions = struct( ...
                         'CaseSensitivity',false, ...
                         'PartialMatching','unique', ...
                         'StructExpand',false, ...
                         'IgnoreNulls',false);
        pstruct  =  coder.internal.parseParameterInputs(params,poptions,args{:});
        ampOut   =  coder.internal.getParameterValue(pstruct.Amplitude,defaultAmp,args{:});
        scaleOut =  coder.internal.getParameterValue(pstruct.Scale,defaultScale,args{:});
        ovpOut   =  coder.internal.getParameterValue(pstruct.OverlapPercent,...
                                             defaultOverlapPercent,args{:});
        winOut   =  coder.internal.getParameterValue(pstruct.Window,defaultWin,args{:});    
    end
    % Do validation of name-value pairs
    amplitude = validatestring(ampOut,{'rms','peak','power'},funName,'Amplitude');
    scale     = validatestring(scaleOut,{'linear','dB'},funName,'Scale');
    validateattributes(ovpOut,{'numeric'},{'scalar','nonnegative','<=',100},...
                       funName,'OP');
    overlapPercent = double(ovpOut(1));
    winList = {'hann','hamming','flattopwin','kaiser','rectwin','chebwin'};
    if iscell(winOut)
        lenWinOut = numel(winOut);
        % If value of 'Window' is a cell, then either it should be a one
        % element cell containing window name or a 2 element cell
        % containing window name and a parameter value
        coder.internal.assert(lenWinOut==1 || lenWinOut == 2,...
                        'signal:rpmmap:InvalidCellForWin');
        win =  validatestring(winOut{1},winList,funName,'Window');
        if lenWinOut == 1 
            % only window name is specified in the cell, set window
            % parameter to [];
            winparam = [];
        else % window paramter is specified.
            % only kaiser and chebwin windows are allowed to have
            % parameters. Assert if parameter is specified for other
            % windows
            coder.internal.assert(strcmpi(win,'kaiser') || ...
              strcmpi(win,'chebwin'),'signal:rpmmap:WindowCellArgumentInvalid',...
              win);           
            if strcmpi(win,'kaiser')
                validateattributes(winOut{2},{'numeric'},...
                   {'real', 'positive', 'finite','scalar'},...
                    funName,'Kaiser beta parameter'); 
            else % The window is chebwin
                 validateattributes(winOut{2},{'numeric'},...
                  {'real', 'positive', 'finite','scalar','>=',45},...
                  funName,'chebwin sidelobe attenuation parameter');
            end
             winparam = winOut{2}(1);
        end
    else
        % If value of window is not a cell, then it should be the name of a
        % valid window. 
         win = validatestring(winOut,winList,funName,'Window');
         winparam = [];
    end
else
    % No name-value pairs are specified. Assign defaults.
    win            =  defaultWin;
    winparam       =  [];
    amplitude      =  defaultAmp;
    scale          =  defaultScale;
    overlapPercent =  defaultOverlapPercent;    
end
end

%--------------------------------------------------------------------------
function [isResolValid,winLength] = validateResolutionValue(resolution,minLength,maxLength,fs,winName,winParam,warn)
%check that resolution falls between minimum and maximum allowable values
%either return a valid resolution or a window length

% Compute resolution of the window
minResValue = enbw(getWindow(winName,maxLength,winParam))*fs/maxLength;
maxResValue = enbw(getWindow(winName,minLength,winParam))*fs/minLength;
sminres = sprintf('%0.5g',minResValue);
smaxres = sprintf('%0.5g',maxResValue);

if resolution <= minResValue
  winLength = maxLength;
  isResolValid = false;
  if(warn)
    coder.internal.warning('signal:rpmmap:ResolutionMustBeInRangeUnder',...
      sminres,smaxres);
  end
elseif resolution >= maxResValue
  winLength = minLength;
  isResolValid = false;
  if(warn)
    coder.internal.warning('signal:rpmmap:ResolutionMustBeInRangeOver',...
      sminres,smaxres);
  end
else
  isResolValid = true;  
  winLength = -1; %Recompute winLength using resolution
end
end

%--------------------------------------------------------------------------
function map = mapAmplitudeScale(mapUnscaled,amplitude,scale,nfft)

% Map amplitude
map = coder.nullcopy(real(mapUnscaled));
switch amplitude
  case 'peak'
    map = oneSidedSpectrum(abs(mapUnscaled),nfft);
    % Scale all components except dc
    map(2:end,:) = sqrt(2)*map(2:end,:);
  case 'rms'
    map = oneSidedSpectrum(abs(mapUnscaled),nfft);
  case 'power'
    map = oneSidedSpectrum(abs(mapUnscaled),nfft).^2;
end

% Map scale
switch scale
  case 'linear'
    %no change is needed
  case 'dB'
    if isequal(amplitude,'power')
      map = 10*log10(map);
    else
      map = 20*log10(map);
    end
end
end

%--------------------------------------------------------------------------
function mp = oneSidedSpectrum(map,nfft)
% Get the correctly scaled one sided spectrum - map has magnitude values
% not power values so scaling should be sqrt(2) instead of 2 (note that we
% already have half the spectrum but we need to scale it correctly).

if isodd(nfft) %fft is odd length
  % Don't double DC component
  mp = [ map( 1, : ); sqrt(2)*map( 2:end, : ) ];
else
  % Don't double DC component or unique Nyquist point
  mp = [ map( 1, : ); sqrt(2)*map( 2:end - 1, : ); map( end, : ) ];
end
end

%--------------------------------------------------------------------------
function win = getWindow(winName,winLength,winParam)
 

switch winName
    case 'kaiser'
        win =  kaiser(winLength,winParam);
    case 'chebwin'
        win = chebwin(winLength,winParam);
    case 'hann'
        win = hann(winLength);
    case 'hamming'
        win = hamming(winLength);
    case 'flattopwin'
        win = flattopwin(winLength);
    otherwise % it is rectwin
        win = rectwin(winLength);
end
end


% LocalWords:  xp resampled fsp fs DFT Omax ordermap nwin noverlap extrap fmax
% LocalWords:  Upsampling rpmordermap rpmfreqmap rms nonsparse paramter
% LocalWords:  sidelobe
