function [aout,eout] = arparest( x, p, method)
%ARPAREST   AR parameter estimation via a specified method.
%   A = ARPAREST(X,ORDER,METHOD) returns the polynomial A corresponding to 
%   the AR parametric signal model estimate of vector X using the specified
%   METHOD.  ORDER is the model order of the AR system.
%
%   Supported methods are: 'covariance' and 'modified' although all of the
%   methods of CORRMTX will work. In particular if 'autocorrelation' is
%   used, the results should be the same as those of ARYULE (but slower).
%
%   [A,E] = ARPAREST(...) returns the variance estimate E of the white noise
%   input to the AR model.

%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        S. Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATION,
%              Prentice-Hall, 1987, Chapter 8.
%        P. Stoica and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
%              Prentice-Hall, 1997, Chapter 3

%   Author(s): R. Losada and P. Pacheco
%   Copyright 1988-2020 The MathWorks, Inc.
%#codegen

narginchk(3,3)
% enforce backwards compatibility with row vectors
xm = signal.internal.toColIfVect(x);
validateattributes(x,{'single','double'},{'nonempty','finite','2d','nonsparse'},...
                      'arparest','X');
validateattributes(p,{'numeric'},{'real','nonnegative','integer','scalar'},...
                      'arparest','ORDER');
p = double(p(1));
coder.internal.assert(strcmp(method,'covariance') || strcmp(method,'modified'),...
                      'signal:arparest:UnknMethod');

% Set up necessary but not sufficient conditions for the correlation
% matrix to be nonsingular. From (Marple)
if coder.const(strcmp(method,'covariance'))
    minlength_x = 2*p;
    multiplier = '3/2';
else % 'modified'
    minlength_x = 3*p/2;
    multiplier = '2';
end
% Do some data sanity testing
if isvector(xm)
    coder.internal.errorIf(size(xm,1) < minlength_x,...
        'signal:arparest:VectorTooSmallForModel','X',multiplier);
else
    coder.internal.errorIf(size(xm,1) < minlength_x,...
        'signal:arparest:MatrixTooSmallForModel','X',multiplier);
end

% 'like' will also copy over the complexity we only want the class
a    = coder.nullcopy(zeros(p+1,size(xm,2),'like',xm));
eout = coder.nullcopy(zeros(1,size(xm,2),class(xm)));

for chan=1:size(xm,2)
   % Generate the appropriate data matrix
   XM = corrmtx(xm(:,chan),p,method);
   Xc = XM(:,2:end);
   X1 = XM(:,1);

   % Coefficients estimated via the covariance method
   a(:,chan) = [1; -Xc\X1];

   % Estimate the input white noise variance
   Cz = X1'*Xc;
   variance = X1'*X1 + Cz*a(2:end,chan);

   % Ignore the possible imaginary part due to numerical errors and force
   % the variance estimate of the white noise to be positive
   eout(:,chan) = abs(real(variance));  
end

aout = a.'; % By convention all polynomials are row vectors

% [EOF] arparest.m
