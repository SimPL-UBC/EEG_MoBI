function [Pxx,F,RPxx,Fc] = computeperiodogram(x,win,nfft,esttype,Fs,options)
%COMPUTEPERIODOGRAM   Periodogram spectral estimation.
%   This function is used to calculate the Power Spectrum Sxx, and the
%   Cross Power Spectrum Sxy.
%
%   Sxx = COMPUTEPERIODOGRAM(X,WIN,NFFT) where x is a vector returns the
%   Power Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Sxy = COMPUTEPERIODOGRAM({X,Y},WIN,NFFT) returns the Cross Power
%   Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Inputs:
%    X           - Signal vector or a cell array of two elements containing
%                  two signal vectors.
%    WIN         - Window
%    NFFT        - Number of frequency points (FFT) or vector of
%                  frequencies at which periodogram is desired
%    ESTTYPE     - A string indicating the type of window compensation to
%                  be done. The choices are: 
%                  'ms'    - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'power' - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'psd'   - compensate for Power Spectral Density (PSD);
%                            maintain correct area under the PSD curve.
%     REASSIGN   - A logical (boolean) indicating whether or not to perform
%                  frequency reassignment
%
%   Output:
%    Sxx         - Power spectrum [Power] over the whole Nyquist interval. 
%      or
%    Sxy         - Cross power spectrum [Power] over the whole Nyquist
%                  interval.
%
%    F           - (vector) list frequencies analyzed
%    RSxx        - reassigned power spectrum [Power] over Nyquist interval
%                  has same size as Sxx.  Empty when 'reassigned' option
%                  not present.
%    Fc          - center of gravity frequency estimates.  Same size as
%                  Sxx.  Empty when 'reassigned' option not present.
%

%   Copyright 1988-2019 The MathWorks, Inc.
%#codegen

narginchk(5,7);
if nargin<6
  reassign = false;
  range = 'twosided';
else
  reassign = options.reassign;
  range = options.range;
end

% use normalized frequencies when Fs is empty
if isempty(Fs) || isnan(Fs)
    Fs = 2*pi;
end

% Validate inputs and convert row vectors to column vectors.
[x1,~,y,is2sig,win1] = validateinputs(x,win,nfft);

% Window the data
xw = bsxfun(@times,x1,win1);


% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels

[Xx,F] = computeDFT(xw,nfft,Fs);
if reassign
  xtw = bsxfun(@times,x1,signal.internal.spectral.dtwin(win1,Fs));
  [Xxc,~] = computeDFT(xtw,nfft,Fs);
  if coder.target('MATLAB')
      Fc = -imag(Xxc ./ Xx);
      Fc(~isfinite(Fc)) = 0;
  else
      siz = size(Xxc);
      num = numel(Xxc);
      Fc = coder.nullcopy(zeros(num,1));
      for i = 1:num
          val = -imag(Xxc(i)/Xx(i));
          if ~isfinite(Fc)
              Fc(i) = 0;
          else
              Fc(i) = val;
          end
      end
      Fc = reshape(Fc,siz);
  end 

  Fc = bsxfun(@plus,F,Fc);
  
else
    Fc=[];
end

% if two signals are used, we are being called from welch and are not
% performing reassignment.
if is2sig
  yw = bsxfun(@times,y,win1);
else
  yw = zeros(0,class(y));
end 

% Evaluate the window normalization constant.  A 1/N factor has been
% omitted since it will cancel below.
if any(strcmpi(esttype,{'ms','power'}))
  if reassign
    if isscalar(nfft)
      U = nfft*(win1'*win1);
    else
      U = numel(win1)*(win1'*win1);
    end
  else
    % The window is convolved with every power spectrum peak, therefore
    % compensate for the DC value squared to obtain correct peak heights.
    U = sum(win1)^2;
  end
else
    U = win1'*win1;  % compensates for the power of the window.
end

if is2sig
  [Yy,F] = computeDFT(yw,nfft,Fs);
  % We use bsxfun here because Yy can be a single vector or a matrix
  Pxx = bsxfun(@times,Xx,conj(Yy))/U;  % Cross spectrum.
else
  % We do a real cast to make Pxx real in codegen.  
  Pxx = real(Xx.*conj(Xx)/U);                % Auto spectrum.
end


% Perform reassignment
if reassign
  RPxx = reassignPeriodogram(Pxx, F, Fc, nfft, range);
else
  RPxx = cast([],'like',Pxx); 
  Fc = [];
end
end



%--------------------------------------------------------------------------
function [x1,Lx,y,is2sig,win1] = validateinputs(x,win,~)
% Validate the inputs to computexperiodogram and convert row vectors to
% column vectors for backwards compatiblility with R2014a and prior
% releases

% Set defaults and convert to row vectors to columns.
is2sig= false;
win1   = win(:);
Lw    = length(win);

% Determine if one or two signal vectors was specified.
if iscell(x)
    if length(x) > 1
        y1 = x{2};
        validateattributes(y1,{'single','double'}, {'finite','nonnan'},'periodogram','x');
        if isvector(y1)
            y = y1(:);
        else
            y=y1;
        end
        is2sig = true;
    end
    x2 = x{1};
else
    y=[];
    x2 = x;
end

validateattributes(x2,{'single','double'}, {'finite','nonnan'},'periodogram','x');
if isvector(x2)
    x1 = x2(:);
else
    x1=x2;
end

Lx = size(x1,1);

if is2sig
    Ly  = size(y,1);
    if Lx ~= Ly
        coder.internal.error('signal:computeperiodogram:invalidInputSignalLength');
    end
    if size(x1,2)~=1 && size(y,2)~=1 && size(x1,2) ~= size(y,2)
        coder.internal.error('signal:computeperiodogram:MismatchedNumberOfChannels');
    end
end

coder.internal.errorIf(Lx ~= Lw, 'signal:computeperiodogram:invalidWindow', 'window');

if (numel(x1)<2 || numel(size(x1))>2)
    coder.internal.error('signal:computeperiodogram:NDMatrixUnsupported');
end

end

% -------------------------------------------------------------------------
function RP = reassignPeriodogram(P, f, fcorr, nfft, range)

% for each column input of Sxx, reassign the power additively
% independently.

nChan = size(P,2);

nf = numel(f);
fmin = f(1);
fmax = f(end);

% compute the destination row for each spectral estimate
% allow cyclic frequency reassignment only if we have a full spectrum
if isscalar(nfft) && strcmp(range,'twosided')
  rowIdx = 1+mod(round((fcorr(:)-fmin)*(nf-1)/(fmax-fmin)),nf);
else
  rowIdx = 1+round((fcorr(:)-fmin)*(nf-1)/(fmax-fmin));
end

if coder.target('MATLAB')

    % compute the destination column for each spectral estimate
    colIdx = repmat(1:nChan,nf,1);

    % reassign the estimates that fit inside the frequency range
    P = P(:);
    idx = find(rowIdx>=1 & rowIdx<=nf);
    RP = accumarray([rowIdx(idx) colIdx(idx)], P(idx), [nf nChan]);
    
else

    rowIdxLen = length(rowIdx);
    newLen=0;
    for i=1:rowIdxLen
        if rowIdx(i) >=1 && rowIdx(i) <= nf
            newLen=newLen+1;
        end
    end

    siz = size(P,1);
    rowIdx1 = coder.nullcopy(zeros(newLen,1));
    p1 = coder.nullcopy(zeros(newLen,1,'like',P));
    channels = coder.nullcopy(zeros(nChan,1));
    tempIdx =1;
    cidx = 1; 
    pidx =1;
    for i=1:rowIdxLen
        if rowIdx(i) >=1 && rowIdx(i) <= nf
            rowIdx1(tempIdx) = rowIdx(i);
            p1(tempIdx) = P(pidx,cidx);
            tempIdx = tempIdx+1;
        end
        pidx = pidx+1;
        if mod(i,siz) == 0
            channels(cidx) = tempIdx -1;
            cidx = cidx+1;
            pidx =1;
        end
    end

    RP = complex(zeros(nf,nChan,'like',P));
    chan =1; 
    for i=1:newLen
         if i>channels(chan)
            chan = chan+1;
        end
        RP(rowIdx1(i),chan) = RP(rowIdx1(i),chan) + p1(i);
    end
end

end

% -------------------------------------------------------------------------

% LocalWords:  Sxx Sxy NFFT ESTTYPE RSxx Fc Fs Yy computexperiodogram
% LocalWords:  compatiblility nonnan
