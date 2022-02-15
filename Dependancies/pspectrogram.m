function varargout = pspectrogram(x,spectype,varargin)
%PSPECTROGRAM Spectrogram and cross spectrogram
%   S = PSPECTROGRAM(X,'spect',WINDOW,NOVERLAP,NFFT,Fs,...)
%   S = PSPECTROGRAM({X,Y},'xspect',WINDOW,NOVERLAP,NFFT,Fs,...)
%   [S,F,T,P,fcorr,tcorr] = PSPECTROGRAM(...)
%
%   Inputs:
%      see "help spectrogram" for complete description of all input
%      arguments. SPECTTYPE is a string specifying the type of estimate to
%      return, the choices are: 'spect' and 'xspect'.
%
%   Outputs:
%      see "help spectrogram" and "help xspectrogram" for complete
%      description of all output arguments.
%      The output definition depends on the input string ESTTYPE:
%      S - complex STFT ('spect') or power cross spectrogram ('xspect')
%      F - frequencies for S
%      T - times for S
%      P - power spectrogram ('spect') or complex cpsd ('xspect')
%      fcorr - matrix of corrected (reassigned) frequencies
%      tcorr - matrix of corrected (reassigned) times

%   Copyright 2016-2020 The MathWorks, Inc.
%#codegen

% Parse arguments specific to PSPECTROGRAM and remove from argument list.
isInMATLAB = coder.target('MATLAB');
[esttype,reassign,faxisloc,threshold,timeDimension,inpArgs] = signal.internal.getSpectrogramOptions(varargin{:});
% Parse input arguments (using the PWELCH parser since it shares the same API).
[xw,nx,~,yw,ny,win,~,~,noverlap,~,~,options] = signal.internal.spectral.welchparse(x,esttype,inpArgs{:});

% Check for valid input signals
chkinput(xw,'X');

% cast to enforce precision rules
noverlap = signal.internal.sigcasttofloat(noverlap,'double',...
    'spectrogram','NOVERLAP','allownumeric');
coder.internal.errorIf(noverlap(1,1)~=round(noverlap(1,1))||noverlap(1,1)<0,...
    'signal:spectrogram:NoverlapNonnegativeInteger');
validateattributes(win,{'double','single'},{'nonsparse'},'spectrogram','Window');

% Process frequency-specific arguments
[Fs,nfft,isnormfreq,options] = processFrequencyOptions(options,reassign);

% Make x and win into columns
xCol = reshape(xw,[],1);
winCol = reshape(win,[],1);
nwin = length(winCol);

% Place x into columns and return the corresponding central time estimates.
[xin,t] = signal.internal.stft.getSTFTColumns(xCol,nx,nwin,noverlap,Fs);

if strcmpi(spectype,'spect') % Spectrogram
    % There is no second input signal, so no need to compute STFT columns.
    yin = [];
    
    % Compute the raw STFT
    % Apply the window to the array of offset signal segments.
    [y1,f] = computeDFT(bsxfun(@times,winCol,xin),nfft,Fs);
    
    % Compute reassignment matrices.
    [fcorr,tcorr] = computeReassign(xin,winCol,Fs,nfft,y1,f,t,reassign,nargout);
 
    
    % Truncate output and adjust any time-frequency corrections based on
    % FREQRANGE spectrum format ('centered', 'onesided', 'twosided')
    [S,f,fcorr,tcorr] = formatSpectrogram(y1,f,fcorr,tcorr,Fs,nfft,options);
    
    % Compute PSD when required.
    if nargout==0 || nargout>3
        [P,f] = compute_PSD(winCol,S,nfft,f,t,options,Fs,esttype,threshold,reassign,fcorr,tcorr);
    else
        P = [];
    end
    
    % Shift the outputs when 'centered' is specified.
    if options.centerdc && length(options.nfft)==1
        [S,f,P,fcorr,tcorr] = centerOutputs(nargout,S,f,P,fcorr,tcorr);
    end
    
    if nargout==0
        coder.internal.assert(isInMATLAB,'signal:spectrogram:plottingNotSupported');
        % plot when no output arguments are specified
        displayspectrogram(t,f,P,isnormfreq,faxisloc,esttype,threshold);
    else
        % Non-conjugate transpose S, t, P, fcorr, and tcorr if
        % timeDimension is 'downrows'
        [St,tt,Pt,fcorrt,tcorrt] = nonConjugateTransposeOutputs(timeDimension,S,t,P,fcorr,tcorr);
    end
    
else % Cross Spectrogram
    % Validate second input argument
    chkinput(yw,'Y');
    yCol = reshape(yw,numel(yw),1);
    
    % Place y into columns.
    yin = signal.internal.stft.getSTFTColumns(yCol,ny,nwin,noverlap,Fs);
    
    if options.centerdc
        freqrange = 'centered';
    else
        freqrange = options.range;
    end
    
    % Compute cross STFT using CPSD.
    [S1,f] = cpsd(xin,yin,winCol,0,nfft,Fs,freqrange);
    % In codegen, cpsd returns a 3D matrix as output if the inputs xin and yin
    % are matrices. Here xin and yin are matrices of the same size, so the
    % output S1 in MATLAB is a 2D matrix. We do the below assignment to
    % enforce the output to be 2D in codegen as well.
    S = S1(:,:,1);
    
    % Scale the cross spectrogram to power using the effective noise
    % bandwidth of the window, if requested.
    if strcmp(esttype,'power')
        S = S*enbw(winCol,Fs);
    end
    
    % Compute time-varying cross spectrum when requested.
    if nargout>3
        P = S;
    else
        P = [];
    end
    
    % Return the cross spectrogram (magnitude) as the first argument.
    S = abs(S);
    
    % Apply threshold
    if threshold>0
        S(S<threshold) = 0;
    end

    if nargout==0
        coder.internal.assert(isInMATLAB,'signal:spectrogram:plottingNotSupported');
        % plot when no output arguments are specified
        displayspectrogram(t,f,S,isnormfreq,faxisloc,esttype,threshold);
    else
        % Non-conjugate transpose S, t, and P if timeDimension is
        % 'downrows'
        [St,tt,Pt] = nonConjugateTransposeOutputs(timeDimension,S,t,P);
    end
    fcorrt = [];
    tcorrt = [];        
end

% cast to enforce precision rules
if nargout > 0
    isSingle =  isa(xin,'single') || isa(yin,'single') || isa(winCol,'single');
    [varargout{1:nargout}] = assignAndcastToSingle(nargout,isSingle,St,f,tt,Pt,fcorrt,tcorrt);
end


%--------------------------------------------------------------------------
function chkinput(x,X)
% Check for valid input signal
coder.internal.errorIf(isempty(x) || issparse(x) || ~isfloat(x),'signal:spectrogram:MustBeFloat', X);
coder.internal.errorIf(~isvector(x),'signal:spectrogram:MustBeVector', X);

%--------------------------------------------------------------------------
function displayspectrogram(t,f,Pxx,isFsnormalized,faxisloc,esttype, threshold)
% Cell array of the standard frequency units strings

if strcmpi(esttype,'power')
    plotOpts.cblbl = getString(message('signal:dspdata:dspdata:PowerdB'));
else
    if isFsnormalized
        plotOpts.cblbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBradsample'));
    else
        plotOpts.cblbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
    end
end

%Threshold in dB
plotOpts.freqlocation = faxisloc;
plotOpts.threshold = 10*log10(threshold+eps);
plotOpts.isFsnormalized = logical(isFsnormalized);

%Power in dB
signalwavelet.internal.convenienceplot.plotTFR(t,f, 10*log10(abs(Pxx)+eps),plotOpts);    

% -------------------------------------------------------------------------
function [Pxx,f,fcorr,tcorr] = compute_PSD(win,y,nfft,f,t,options,Fs,esttype,threshold,reassign,fcorr,tcorr)

% Evaluate the window normalization constant.
if strcmpi(esttype,'power')
    if reassign
        % compensate for the power of the window including a
        % 1/N scaling factor omitted by FFT/DFT computation.
        if isscalar(nfft)
            U = nfft(1)*(win'*win);
        else
            U = numel(win)*(win'*win);
        end
    else
        % The window is convolved with every power spectrum peak, therefore
        % compensate for the DC value squared to obtain correct peak heights.
        % The 1/N factor has been omitted since it will cancel below.
        U = sum(win)^2;
    end
else
    % compensates for the power of the window.
    % The 1/N factor has been omitted since it will cancel below.
    U = win'*win;
end

Sxx = real(y.*conj(y))/U; % Auto spectrum.

% reassign in-place when requested.
if reassign
    Sxx1 = signal.internal.spectral.reassignSpectrum(Sxx, f, t, fcorr, tcorr, options);
else
    Sxx1 = Sxx;
end

% Compute the one-sided or two-sided PSD [Power/freq]. Also compute
% the corresponding half or whole power spectrum [Power].
[Pxx,f] = computepsd(Sxx1, f, options.range, nfft, Fs, esttype);

% remove low-power estimates if requested
if threshold>0
    Pxx(Pxx<threshold) = 0;
end

% -------------------------------------------------------------------------
function [y,f,Pxx,fcorr,tcorr] = centerOutputs(nOut,y,f,Pxx,fcorr,tcorr)
% center y,fcorr,tcorr only if specified in the output list
% center f,Pxx if specified in the output list or when plotting
% nOut contains the number of output arguments of SPECTROGRAM
if nOut>0
    y = signal.internal.spectral.centerest(y);
end

if nOut==0 || nOut>1
    f = signal.internal.spectral.centerfreq(f);
end

if nOut==0 || nOut>3
    Pxx = signal.internal.spectral.centerest(Pxx);
end

if nOut>4
    fcorr = signal.internal.spectral.centerest(fcorr);
end

if nOut>5
    tcorr = signal.internal.spectral.centerest(tcorr);
end

% -------------------------------------------------------------------------
function varargout = assignAndcastToSingle(nOut,isSingle,varargin)
% convert outputs to single precision when specified in output list
% nOut contains the number of output arguments of SPECTROGRAM
out = cell(1,nOut);
for i = 1:nOut
    if isSingle
        out{i} = single(varargin{i});
    else
        out{i} = varargin{i};
    end  
end
[varargout{1:nOut}] = out{:};  

% -------------------------------------------------------------------------
function varargout = nonConjugateTransposeOutputs(timeDimension,varargin)
% Non-conjugate transpose outputs if timeDimension is not the default
% value.

out = cell(1,nargout);
for i=1:nargin-1
    if strcmpi(timeDimension,'acrosscolumns')
        out{i} = varargin{i};
    else
        out{i} = varargin{i}.';
    end
end
[varargout{1:nargout}] = out{:};

% -------------------------------------------------------------------------
function [Fs,nfft,isnormfreq,options] = processFrequencyOptions(options,reassign)
% Determine whether an empty was specified for Fs (i.e., Fs=1Hz) or
% returned by welchparse which means normalized Fs is used.

% Cast to enforce Precision rules
Fs = double(options.Fs(1));
nfft = double(options.nfft);

% when Fs is specified as [], welchparse() returns 1 Hz.
% welchparse() returns nan only when Fs is omitted
isnormfreq = isnan(Fs);
if isnormfreq
    Fs = 2*pi;
end

if length(nfft) > 1
    % Frequency vector was specified, return and plot two-sided PSD
    if strcmpi(options.range,'onesided')
        coder.internal.warning('signal:welch:InconsistentRangeOption');
    end
    options.range = 'twosided';
end

% prevent unneeded temporary conversion to one-sided spectrum
if options.centerdc
    options.range = 'twosided';
end

% ensure frequency vector is linearly spaced when performing reassignment
if reassign && ~isscalar(nfft)
    f = options.nfft(:);
    
    % see if we can get a uniform spacing of the freq vector
    [~, ~, ~, maxerr] = signal.internal.spectral.getUniformApprox(f);
    
    % see if the ratio of the maximum absolute deviation relative to the
    % largest absolute in the frequency vector is less than a few eps
    if options.isNFFTSingle
        isuniform = maxerr < 3*eps('single');
    else
        isuniform = maxerr < 3*eps('double');
    end
    if ~isuniform
        coder.internal.error('signal:spectrogram:ReassignFreqMustBeUniform');
    end
end


% -------------------------------------------------------------------------
function [yout,fout,fcorrout,tcorrout]  = formatSpectrogram(y,f,fcorr,tcorr,Fs,nfft,options)
% truncate output and adjust any time-frequency corrections based on
% FREQRANGE spectrum format ('centered', 'onesided', 'twosided')

% if nfft is a scalar, it is the length of the fft, otherwise it contains
% the output frequency vector
freqvecspecified = length(nfft)>1;
% truncate output when using one-sided spectrum
if ~freqvecspecified && strcmpi(options.range,'onesided')
    fout = psdfreqvec('npts',nfft,'Fs',Fs,'Range','half');
    yout = y(1:length(fout),:);
    if ~isempty(fcorr)
        fcorrout = fcorr(1:length(fout),:);
    else
        fcorrout = fcorr;
    end
    if ~isempty(tcorr)
        tcorrout = tcorr(1:length(fout),:);
    else
        tcorrout = tcorr;
    end
else
    yout = y;
    fcorrout = fcorr;
    tcorrout = tcorr;
    fout = f;
end

if ~isempty(fcorrout)
    if options.centerdc || freqvecspecified && any(nfft < 0,'all')
        % map to [-Fs/2,Fs/2) when using negative frequencies
        fcorrout = mod(fcorrout+Fs/2,Fs)-Fs/2;
    else
        % map to [0,Fs) when using positive frequencies
        fcorrout = mod(fcorrout,Fs);
    end
end

%--------------------------------------------------------------------------
function [fcorr,tcorr] = computeReassign(xin,win,Fs,nfft,y,f,t,reassign,narg)
% Compute the reassignment time and frequency matrices fcorr and tcorr
if reassign || narg>4
    % Apply frequency correction from time derivative window
    yc = computeDFT(bsxfun(@times,signal.internal.spectral.dtwin(win,Fs),xin),nfft,Fs);
    fcorr = -imag(yc ./ y);
    fcorr(~isfinite(fcorr)) = 0;
    fcorr = bsxfun(@plus,f,fcorr);
else
    fcorr = [];
end

if reassign || narg>5
    % Apply time correction from frequency derivative window
    yc = computeDFT(bsxfun(@times,signal.internal.spectral.dfwin(win,Fs),xin),nfft,Fs);
    tcorr = real(yc ./ y);
    tcorr(~isfinite(tcorr)) = 0;
    tcorr = bsxfun(@plus,t,tcorr);
else
    tcorr = [];
end

% LocalWords:  spect NOVERLAP NFFT Fs xspect fcorr tcorr SPECTTYPE xspectrogram
% LocalWords:  ESTTYPE STFT allownumeric welchparse Noverlap FREQRANGE downrows
% LocalWords:  xin Powerd Powerfrequencyd Bradsample BHz DFT Pxx acrosscolumns
% LocalWords:  nfft npts
