function [Pxx, F, Frange, rbw, extraArgs, status] = psdparserange(funcName, beta, varargin)
%PSDPARSERANGE Helper function for ranged power and psd estimates
%   computes and/or returns a PSD for use based upon the input
%
%      funcName - name of calling function to use with error messages
%
%      beta - coefficient for Kaiser window.
%             Usually 0 (rectangular) or 38 (for use with sinusoids)
%   
%      varargin - supports:
%                    FUNC(X) 
%                    FUNC(X, Fs)
%                    FUNC(Pxx, F)
%                    FUNC(Sxx, F, RBW)
%                    FUNC(..., F(or Fs), FREQRANGE)
%                    FUNC(..., F(or Fs), FREQRANGE, EXTRAARGS)
%          X  - time-domain input
%          Fs - sample rate
%          F - frequency vector
%          RBW - resolution bandwidth
%          Pxx - PSD estimate
%          Sxx - spectral power estimate
%          FREQRANGE - must be empty or a two-element vector
%          EXTRAARGS - anything following the FREQRANGE
%          NORMF - normalized frequency flag

%   Copyright 2014-2019 The MathWorks, Inc.
%#codegen

n = numel(varargin);

oneSided = false;
hasNyquist = false;

if n==1 || ...
        (coder.internal.isConst(isempty(varargin{2}) || isscalar(varargin{2})) ...
        && (isempty(varargin{2}) || isscalar(varargin{2})))
    % if varargin{2} is a scalar or empty, it represents Fs, else F
    inputType = 'time';
    oneSided = true;
    [Pxx, F, Frange, rbw] = parseTime(funcName, beta, varargin{1:min(n,3)});
    hasNyquist = hasNyquistBin(varargin{1});
    extraArgs = {varargin{4:end}};
elseif n==2
    inputType = 'psd';
    [Pxx, F, Frange, rbw] = parsePSD(funcName, varargin{1:min(n,3)});
    extraArgs = {varargin{4:end}};
else
    if (coder.internal.isConst(isscalar(varargin{3})) ...
            && (isscalar(varargin{3})))
        % if varargin{3} is a scalar, it represents RBW, else FREQRANGE
        inputType = 'power';
        [Pxx, F, Frange, rbw] = parsePower(funcName, varargin{1:min(n,4)});
        extraArgs = {varargin{5:end}};
    else
        inputType = 'psd';
        [Pxx, F, Frange, rbw] = parsePSD(funcName, varargin{1:min(n,3)});
        extraArgs = {varargin{4:end}};
    end
end

% report normalized frequency flag
status.normF = (n==1) || isempty(varargin{2});
status.inputType = inputType;
status.oneSided = oneSided;
status.hasNyquist = hasNyquist;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function flag = hasNyquistBin(x)
if isvector(x)
  flag = iseven(numel(x));
else
  flag = iseven(size(x,1));
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [Pxx, F, Frange, rbw] = parseTime(funcName, beta, x, fs, Frange)

% validate x
validateattributes(x,{'numeric'},{'2d','finite'},funcName,'X',1);
if isrow(x)
    x_t = x.';
else
    x_t = x;
end

% validate fs
if nargin < 4
    % use angular frequency
    fs_t = 2*pi;
else
    if isempty(fs)
        % use angular frequency
        fs_t = cast(2*pi,'like',fs);
    else
        validateattributes(fs,{'numeric'},{'real','finite','scalar','positive'}, ...
            funcName,'Fs',2);
        % scalar inference in codegen
        fs_t = fs(1);
    end
end

% use Kaiser window to reduce effects of leakage
n = size(x_t,1);
w = kaiser(n,beta);
rbw = enbw(w,fs_t);

% use one-sided PSD for real signals, otherwise use centered for complex
if isreal(x_t)
    [Pxx,F] = periodogram(x_t,w,n,fs_t,'psd');
else
    [Pxx,F] = periodogram(x_t,w,n,fs_t,'centered','psd');
end

% check freq range vector. If not specified, return empty.
if nargin < 5
    Frange = [];
else
    if ~isempty(Frange)
        if fs_t==2*pi
            validateattributes(Frange,{'numeric'},{'real','finite','increasing','size',[1 2]}, ...
                funcName,'FREQRANGE',2);
        else
            validateattributes(Frange,{'numeric'},{'real','finite','increasing','size',[1 2]}, ...
                funcName,'FREQRANGE',3);
        end
        coder.internal.errorIf(Frange(1)<F(1) || Frange(2)>F(end),'signal:psdparserange:FreqRangeOutOfBounds');
    end
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [PxxOut, FOut, Frange, rbw] = parsePSD(funcName, Pxx, F, Frange)

% validate Pxx
validateattributes(Pxx,{'numeric'},{'2d','real','nonnegative'},funcName,'Pxx',1);
if isrow(Pxx)
    PxxOut = Pxx.';
else
    PxxOut = Pxx;
end

% validate F
validateattributes(F,{'numeric'},{'real','vector','finite'},funcName,'F',2);
% F should not become a scalar at runtime if code was generated using
% var-sized inputs as a scalar represents Fs
coder.internal.errorIf(isscalar(F),'signal:psdparserange:FCannotBecomeScalar');
if isrow(F)
    FOut = F.';
else
    FOut = F;
end
coder.internal.assert(size(PxxOut,1) == numel(FOut),'signal:psdparserange:FreqVectorMismatch');

% check freq range vector. If not specified, return empty.
if nargin < 4
    Frange = [];
else
    % Frange should not become a scalar at runtime if code was generated
    % using var-sized inputs as a scalar represents RBW
    coder.internal.errorIf(isscalar(Frange),'signal:psdparserange:FreqRangeCannotBecomeScalar');
    if ~isempty(Frange)
        validateattributes(Frange,{'numeric'},{'real','finite','increasing','size',[1 2]}, ...
            funcName,'FREQRANGE',3);
        coder.internal.errorIf(Frange(1)<FOut(1) || Frange(2)>FOut(end),'signal:psdparserange:FreqRangeOutOfBounds');
    end
end

% no need to return an RBW for PSD
rbw = NaN;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [Pxx, FOut, Frange, rbwOut] = parsePower(funcName, Sxx, F, rbw, Frange)

% validate Sxx
validateattributes(Sxx,{'numeric'},{'2d','real','nonnegative'},funcName,'Sxx',1);
if isrow(Sxx)
    Sxx_t = Sxx.';
else
    Sxx_t = Sxx;
end

% validate F
validateattributes(F,{'numeric'},{'real','vector','finite'},funcName,'F',2);
if isrow(F)
    FOut = F.';
else
    FOut = F;
end
coder.internal.assert(size(Sxx_t,1) == numel(FOut),'signal:psdparserange:FreqVectorMismatch');

% ensure specified RBW is larger than a bin width
df = mean(diff(FOut,[],1),1);
validateattributes(rbw,{'numeric'},{'real','finite','positive','scalar','>=',df}, ...
    funcName,'RBW',3);
rbwOut = rbw(1);

% check freq range vector. If not specified, return empty.
if nargin < 5
    Frange = [];
else
    if ~isempty(Frange)
        validateattributes(Frange,{'numeric'},{'real','finite','increasing','size',[1 2]}, ...
            funcName,'FREQRANGE',3);
        coder.internal.errorIf(Frange(1)<FOut(1) || Frange(2)>FOut(end),'signal:psdparserange:FreqRangeOutOfBounds');
    end
end

Pxx = Sxx_t./rbwOut;

end