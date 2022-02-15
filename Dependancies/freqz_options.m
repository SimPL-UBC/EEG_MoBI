function [nfft,fvflag,range,Fs] = freqz_options(varargin)
%#codegen

%   Copyright 2020 The MathWorks, Inc.

narginchk(0,3);
if nargin == 0
    nfft = 512;
    fvflag = 0;
    range = 'onesided';
    Fs = [];
else
    nHalf = 0;
    nWhole = 0;
    nOnesided = 0;
    nTwosided = 0;
    idxChar = 0;
    for i = 1:length(varargin)
        coder.internal.assert(isnumeric(varargin{i}) || ischar(varargin{i}) ||....
            isStringScalar(varargin{i}),'signal:signalanalysisbase:invalidInput');
        if ischar(varargin{i}) || isStringScalar(varargin{i})
            coder.internal.assert(coder.internal.isConst(varargin{i}),...
                'signal:signalanalysisbase:inputStringNotConstant');
            coder.internal.assert(~isempty(varargin{i}) && ...
                (strncmpi(varargin{i},'whole',strlength(varargin{i})) ||...
                 strncmpi(varargin{i},'half',strlength(varargin{i})) || ...
                 strncmpi(varargin{i},'onesided',length(varargin{i})) || ...
                 strncmpi(varargin{i},'twosided',length(varargin{i}))),...
                 'signal:signalanalysisbase:invalidString');
            
            if coder.const(strncmpi(varargin{i},'half',strlength(varargin{i})))
                nHalf = nHalf + 1;
            elseif coder.const(strncmpi(varargin{i},'whole',strlength(varargin{i})))
                nWhole = nWhole + 1;
            elseif coder.const(strncmpi(varargin{i},'onesided',strlength(varargin{i})))
                nOnesided = nOnesided + 1;
            elseif coder.const(strncmpi(varargin{i},'twosided',strlength(varargin{i})))
                nTwosided = nTwosided + 1;
            end           
            idxChar = i;
        end
    end
    
    coder.internal.assert(nHalf + nWhole + nOnesided + nTwosided <=1,...
                'signal:signalanalysisbase:duplicateStrings');
    if coder.const(nWhole || nTwosided)
        range = 'twosided';
    else
        % either 'onesided'/'half' was specified or no string inputs were
        % given.
        range = 'onesided';
    end
    
    if coder.const(idxChar > 0)
        args = {varargin{1:idxChar-1},varargin{idxChar+1:end}};
    else
        args = varargin;
    end
    
    if ~isempty(args)
        if isempty(args{1})
            nfft = 512;
        else
            validateattributes(args{1},{'numeric'},...
                {'real','vector'},'','n');
            nfft = double(args{1});
        end
        if length(args) > 1
            if isempty(args{2})
                Fs = 1;
            else
                validateattributes(args{2},{'numeric'},...
                    {'scalar','real','positive'},'','Fs');
                Fs = double(args{2}(1));
            end
        else
            Fs = [];
        end
    else
        nfft = 512;
        Fs = [];
    end
    
    isNFFT = coder.internal.isConst(isscalar(nfft)) && isscalar(nfft);
    if isNFFT
        fvflag = 0;
    else
        if isscalar(nfft)
            coder.internal.error(...
                'signal:signalanalysisbase:varSizeNfftCannotBecomeScalar')
        end
        fvflag = 1;
    end
end

% LocalWords:  signalanalysisbase Fs Nfft
