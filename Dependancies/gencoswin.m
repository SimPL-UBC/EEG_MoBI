function w = gencoswin(window,N,samplingFlag)
%GENCOSWIN   Returns one of the generalized cosine windows.
%   GENCOSWIN returns the generalized cosine window specified by the
%   first string argument. Its inputs can be
%     Window name    - a string, any of 'hamming', 'hann', 'blackman'.
%     N              - length of the window desired.
%     Sampling flag  - optional string, one of 'symmetric', 'periodic'.

%   Copyright 1988-2018 The MathWorks, Inc.

%#codegen

% number of arguments check
narginchk(2,3);

% Cast to enforce Precision Rules
validateattributes(N,{'numeric'},{'real','finite'},window,'N',1);
N = double(N);

% Check for trivial orders
[L,w,trivialwin] = check_order(N);
if trivialwin 
    return
end

% Select the sampling option
if nargin == 2 % no sampling flag specified, use default.
    sflag = 'symmetric';
else
    validateattributes(samplingFlag,{'char','string'},{'scalartext'});
    sflag = convertStringsToChars(samplingFlag);
end

% Allow partial strings for sampling options
allsflags = {'symmetric','periodic'};
sflag = validatestring(sflag,allsflags,window,'sflag',2);

% Evaluate the window
switch sflag
    case 'periodic'
        w = sym_window(L+1,window,2);        
    case 'symmetric'
        w = sym_window(L,window,1);
end

%---------------------------------------------------------------------
function w = sym_window(n,window,startIdx)
%SYM_WINDOW   Symmetric generalized cosine window.
%   SYM_WINDOW Returns an exactly symmetric N point generalized cosine
%   window by evaluating the first half and then flipping the same samples
%   over the other half.

if iseven(n)
    % Even length window
    half = n/2;
    w = calc_window(half,n,window);
    w = [w; w(end:-1:startIdx)];
else
    % Odd length window
    half = (n+1)/2;
    w = calc_window(half,n,window);
    w = [w; w(end-1:-1:startIdx)];
end

%---------------------------------------------------------------------
function w = calc_window(m,n,window)
%CALC_WINDOW   Calculate the generalized cosine window samples.
%   CALC_WINDOW Calculates and returns the first M points of an N point
%   generalized cosine window determined by the 'window' string.

x = (0:m-1)'/(n-1);

switch window
    case 'hann'
        % Hann window
        % w = 0.5 * (1 - cos(2*pi*(0:m-1)'/(n-1)));
        w = 0.5 - 0.5*cos(2*pi*x);
    case 'hamming'
        % Hamming window
        % w = (54 - 46*cos(2*pi*(0:m-1)'/(n-1)))/100;
        w = 0.54 - 0.46*cos(2*pi*x);
    case 'blackman'
        % Blackman window
        % Force end points to zero to avoid close-to-zero negative values caused
        % by roundoff errors.
        % w = (42 - 50*cos(2*pi*(0:m-1)/(n-1)) + 8*cos(4*pi*(0:m-1)/(n-1)))'/100;
        w = 0.42 - 0.5*cos(2*pi*x) + 0.08*cos(4*pi*x);
        if ~isempty(w) % to take care of special case [], this check is required.
            w(1) = 0;
        end
    case 'flattopwin'
        % Flattop window
        % Coefficients as defined in the reference [1] (see flattopwin.m)
        a0 = 0.21557895;
        a1 = 0.41663158;
        a2 = 0.277263158;
        a3 = 0.083578947;
        a4 = 0.006947368;
        w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + ...
            a4*cos(8*pi*x);
end

% [EOF] gencoswin.m

% LocalWords:  allownumeric scalartext sflag CALC roundoff
