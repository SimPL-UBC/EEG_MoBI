function options = musicOptions(x,p,varargin)
%MUSIC_OPTIONS   Parse the optional inputs to the MUSIC function.
%   MUSIC_OPTIONS returns a structure, OPTIONS, with the following fields:
%
%   options.nfft         - number of freq. points at which the psd is estimated
%   options.Fs           - sampling freq. if any
%   options.range        - 'onesided' or 'twosided' pseudospectrum (they correspond to
%                          'half' and 'whole' respectively, but are returned as is by
%                          psdoptions.m
%   options.nw           - number of columns in the data matrix
%   options.noverlap     - number of samples to overlap
%   options.window       - a vector with window coefficients
%   options.CorrFlag     - a flag indicating whether the input is a correlation matrix
%   options.EVFlag       - flag, 0 = MUSIC method ; 1 = EigenVector method
%   options.CorrMatrOrd  - order of the correlation matrix to be used in computations
%#codegen

if coder.target('MATLAB')
    xIsReal = isreal(x);
    % Assign Defaults
    options.nw = [];
    options.noverlap = [];
    options.window = [];
    options.nfft = 256;
    options.Fs = [];
    options.CorrFlag = 0;
    options.EVFlag = 0;
    options.conflevel = 'omitted'; %Default initial value

    % Determine if frequency vector specified
    freqVecSpec = false;
    if (~isempty(varargin) && isnumeric(varargin{1}) && length(varargin{1}) > 1)
        freqVecSpec = true;
    end

    if xIsReal && ~freqVecSpec
        options.range = 'onesided';
    else
        options.range = 'twosided';
    end

    [options,msg,msgobj] = psdoptions(xIsReal,options,varargin{:});

    if isfield(options,'conflevel') && ~strcmp(options.conflevel, 'omitted')
        error(message('signal:music:UnsupportedConfidenceLevels'));
    end

    if length(options.nfft) > 1
        if strcmpi(options.range,'onesided')
            warning(message('signal:music:InconsistentRangeOption'));
        end
        options.range = 'twosided';
    end
    options.CorrMatrOrd = 2*p(1);
    if ~isempty(msg)
        error(msgobj);
    end
else
    coder.internal.prefer_const(p,varargin);
    one  = coder.internal.indexInt(1);
    zero = coder.internal.indexInt(0);
    isrealx = isreal(x);
    validStrings = {'onesided','twosided','centered','ev','corr','whole','half'};
    nOnesided = zero;
    nTwosided = zero;
    nCentered = zero;
    nEv       = zero;
    nCorr     = zero;
    nWhole    = zero;
    nHalf     = zero;
    nNumeric  = zero;
    coder.unroll();
    for i = one:length(varargin)
        coder.internal.assert((ischar(varargin{i}) || isStringScalar(varargin{i})) ||...
            (isnumeric(varargin{i}) && ~issparse(varargin{i}) && ...
            isreal(varargin{i}) && any(size(varargin{i}) <= 1)),...
            'signal:psdoptions:NeedValidOptionsType');
   
        if ischar(varargin{i}) || isStringScalar(varargin{i})
            coder.internal.assert(coder.internal.isConst(varargin{i}),...
                'signal:codegeneration:OptionInputsMustBeConstant')
            coder.const(validatestring(varargin{i},validStrings,'music'));
            if coder.const(strncmpi(varargin{i},'onesided',strlength(varargin{i})))
                nOnesided = nOnesided + 1;
            elseif coder.const(strncmpi(varargin{i},'twosided',strlength(varargin{i})))
                nTwosided = nTwosided + 1;
            elseif coder.const(strncmpi(varargin{i},'whole',strlength(varargin{i})))
                nWhole = nWhole + 1;
            elseif coder.const(strncmpi(varargin{i},'half',strlength(varargin{i})))
                nHalf = nHalf + 1;
            elseif coder.const(strncmpi(varargin{i},'centered',strlength(varargin{i})))
                nCentered = nCentered + 1;
            elseif coder.const(strncmpi(varargin{i},'corr',strlength(varargin{i})))
                nCorr = nCorr + 1;
            elseif coder.const(strcmpi(varargin{i},'ev'))
                nEv = nEv + 1;
            end
        else % we have a numeric option
            nNumeric = nNumeric + one;
            if coder.const(nNumeric == 1)
                if isempty(varargin{i})
                    nfft = 256;
                else
                    nfft = double(varargin{i});
                end
            elseif coder.const(nNumeric == 2)
                if isempty(varargin{i})
                    Fs = 1;
                else
                    validateattributes(varargin{i},{'numeric'},{'scalar',...
                        'positive','finite'},'music','SAMPLE RATE')
                    Fs = double(varargin{i}(1));
                end
            elseif coder.const(nNumeric == 3)
                if coder.internal.isConst(isempty(varargin{i})) && isempty(varargin{i})
                    nw = [];
                    window = [];
                elseif coder.internal.isConst(isscalar(varargin{i})) && isscalar(varargin{i})
                    nw = varargin{i};
                    window = [];
                else
                    % either a fixed size vector of a variable-sized input
                    % at compile time. This will be treated as the window
                    % vector. Error will be thrown if window becomes scalar
                    % or [] at runtime.
                    if isempty(varargin{i}) || isscalar(varargin{i})
                        coder.internal.error(...
                            'signal:codegeneration:VarsizeInputCannotBecomeScalarOrEmpty',...
                            'nw');
                    end
                    window = varargin{i};
                    nw = length(window);
                end
            elseif coder.const(nNumeric == 4)
                if isempty(varargin{i})
                    if isempty(nw)
                        noverlap = [];
                    else
                        noverlap = nw - 1;
                    end
                else
                    validateattributes(varargin{i},{'numeric'},{'scalar',...
                        'finite','nonnegative','integer'},'music','NOVERLAP')
                    noverlap = double(varargin{i}(1));
                end
            end
        end
    end
    coder.internal.assert(nNumeric <= 4,...
        'signal:psdoptions:TooManyNumericOptions');
    if coder.const(nNumeric == 0)
        % No numeric options specified
        Fs       =  [];
        nfft     =  256;
        nw       =  [];
        window   =  [];
        noverlap =  [];
    elseif coder.const(nNumeric == 1)
        % only nfft was specified
        Fs = [];
        nw = [];
        window   =  [];
        noverlap =  [];
    elseif coder.const(nNumeric == 2)
        % nfft and Fs were specified
        nw = [];
        window   =  [];
        noverlap =  [];
    elseif coder.const(nNumeric == 3)
        % nfft, Fs and window/nw were specified
        if isempty(nw)
            noverlap = [];
        else
            noverlap = nw - 1;
        end
    end

    % A given string must be specified only once.
    % Since 'whole' is same as 'twosided', these options should not be specified
    % together. The same holds for 'onesided' and 'half'.
    coder.internal.assert(nOnesided + nHalf  <= 1,'signal:psdoptions:MultipleValues');
    coder.internal.assert(nTwosided + nWhole <= 1,'signal:psdoptions:MultipleValues');
    coder.internal.assert(nCentered  <= 1,'signal:psdoptions:MultipleValues');
    coder.internal.assert(nCorr  <= 1,'signal:psdoptions:MultipleValues');
    coder.internal.assert(nEv  <= 1,'signal:psdoptions:MultipleValues');

    % The following combinations are illegal (conflicting):
    % 1) {'onesided','twosided'},{'half','whole'},{'onesided','whole'},{'half','twosided'}
    % 2) {'onesided','centered'},{'half','centered'}
    coder.internal.assert(nOnesided + nTwosided <= 1,'signal:psdoptions:ConflictingOptions',...
        'onesided','twosided');

    coder.internal.assert(nHalf + nWhole <= 1,'signal:psdoptions:ConflictingOptions',...
        'half','whole');

    coder.internal.assert(nHalf + nTwosided <= 1,'signal:psdoptions:ConflictingOptions',...
        'half','twosided');

    coder.internal.assert(nOnesided + nWhole <= 1,'signal:psdoptions:ConflictingOptions',...
        'onesided','whole');

    coder.internal.assert(nCentered + nOnesided <= 1,'signal:psdoptions:ConflictingOptions',...
        'centered','onesided');

    coder.internal.assert(nCentered + nHalf <= 1,'signal:psdoptions:ConflictingOptions',...
        'centered','half');

    isNFFTScalar = coder.internal.isConst(isscalar(nfft)) && isscalar(nfft);
    if ~isNFFTScalar
        % variable sized nfft is treated as the frequency vector. Error out if
        % it becomes scalar/empty at runtime. We already assigned 256 to nfft
        % if it became empty at runtime. So only a scalar check is required.
        if isscalar(nfft)
            coder.internal.error(...
                'signal:codegeneration:VarsizeInputCannotBecomeScalarOrEmpty','nfft');
        end
    end
 
    if coder.const(nOnesided || nHalf)
        coder.internal.assert(isrealx,'signal:psdoptions:ComplexInputDoesNotHaveOnesidedPSD');
        % ignore the 'onesided' option if a frequency vector is specified
        if ~isNFFTScalar
            coder.internal.warning('signal:music:InconsistentRangeOption');
            range = 'twosided';
        else
            range = 'onesided';
        end
    elseif coder.const(nTwosided || nWhole)
        range = 'twosided';
    else % nothing specified for range.
        if isrealx && isNFFTScalar
            range = 'onesided';
        else
            range = 'twosided';
        end
    end

    if coder.const(nCentered)
        % 'centered' option requires nfft to be a scalar.
        coder.internal.assert(isNFFTScalar,'signal:music:CannotCenterFrequencyVector');
        centerdc  = true;
    else
        centerdc  = false;
    end

    if coder.const(nEv)
        EVFlag = true;
    else
        EVFlag = false;
    end

    if coder.const(nCorr)
        CorrFlag = true;
    else
        CorrFlag = false;
    end

    options = ...
        coder.internal.stickyStruct('CorrMatrOrd',2*p(1),...
        coder.internal.stickyStruct('centerdc',centerdc,...
        coder.internal.stickyStruct('range',range,...
        coder.internal.stickyStruct('EVFlag',EVFlag,...
        coder.internal.stickyStruct('CorrFlag',CorrFlag,...
        coder.internal.stickyStruct('Fs',Fs,...
        coder.internal.stickyStruct('nfft',nfft,...
        coder.internal.stickyStruct('window',window,...
        coder.internal.stickyStruct('noverlap',noverlap,...
        coder.internal.stickyStruct('nw',nw,...
        coder.internal.stickyStruct()))))))))));
end
