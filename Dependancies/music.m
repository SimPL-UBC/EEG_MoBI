function music_data = music(x,p,varargin)
%MUSIC  Implements the heart of the MUSIC algorithm of line spectra estimation.
%   MUSIC is called by both PMUSIC and ROOTMUSIC.
%   
%   Inputs:
%     
%     x - vector or matrix. If vector it is a signal, if matrix it may be either a data
%         matrix such that x'*x=R, or a correlation matrix R.
%     p - scalar or two element vector. If scalar, it indicates the dimension of the 
%         signal subspace. If vector, p(2) is a threshold used to determine the 
%         aforementioned dimension.
%     nfft - (optional) to be used only with PMUSIC. A scalar indicating the number of
%            points used in the evaluation of the pseudospectrum.
%     Fs - (optional) a scalar specifying the sampling frequency. If omitted, we work
%          in rad/sample; if empty it defaults to 1 Hz.
%     nw - (optional) a scalar or vector indicating either the order of the correlation
%          matrix or (when a vector) a window whose length is the order of the matrix
%          and whose values are used to window each column of the data matrix.
%     noverlap - (optional) a integer indicating the number of samples to overlap from
%                column to column.
%     strings - Optional input strings are: 'corr', 'EV' and range ('half' or 'whole').
%
%   Outputs:
%
%     msg - a possible error message.
%
%     music_data - a structure with the following fields:
%          
%          noise_eigenvects - a matrix whose columns are the noise subspace eigenvectors.
%          signal_eigenvects - a matrix whose columns are the signal subspace eigenvectors.
%          eigenvals - the eigenvalues of the correlation matrix.
%          p_eff - the effective dimension of the signal subspace.
%          nfft - number of points used to evaluate the pseudospectrum (only used in PMUSIC).
%          Fs - sampling freq.
%          range - string indicating whether 'half' or the 'whole' pseudospectrum should be
%                  computed. (Only used in PMUSIC.)
%          centerdc - true if 'centered' is specified
%          EVFlag - flag, 0 = MUSIC method; 1 = EigenVector method.

%   Author(s): R. Losada
%   Copyright 1988-2017 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] S. J. Orfanidis, Optimum Signal Processing. An Introduction. 
%         2nd Ed., Macmillan, 1988.
%#codegen
isMATLAB = coder.target('MATLAB');

coder.internal.assert(~isempty(p),'signal:music:EmptySubspaceDimension');
coder.internal.assert(numel(p) <= 2 ,'signal:music:PMustBeVectorOfOneOrTwoElements');
coder.internal.assert(imag(p(1)) == 0 && round(p(1)) == p(1),...
    'signal:music:SubspaceDimensionMustBeInteger');

coder.internal.assert(numel(p) <= 1 || imag(p(2)) == 0,...
    'signal:music:SubspaceThresholdMustBeReal')
p1 = double(real(p));
opts = musicOptions(x,p1,varargin{:});

% Compute the eigenvalues and eigenvectors of the correlation matrix
[eigenvals,eigenvects] = computeeig(x,opts.CorrFlag,opts.CorrMatrOrd,opts.nw,opts.noverlap,opts.window,opts.EVFlag);

% Determine the effective dimension of the signal subspace
p_eff = determine_signal_space(p1,eigenvals);

% Separate the signal and noise eigenvectors
signal_eigenvects = eigenvects(:,1:p_eff);
noise_eigenvects = eigenvects(:,p_eff+1:end);

% Generate the output structure
if isMATLAB
    music_data.noise_eigenvects = noise_eigenvects;
    music_data.signal_eigenvects = signal_eigenvects;
    music_data.eigenvals = eigenvals;
    music_data.p_eff = p_eff;
    music_data.nfft = opts.nfft;
    music_data.Fs = opts.Fs;
    music_data.EVFlag = opts.EVFlag;
    music_data.range = opts.range;
    music_data.centerdc = opts.centerdc;
else
    music_data = coder.internal.stickyStruct('centerdc',opts.centerdc,...
        coder.internal.stickyStruct('range',opts.range,...
        coder.internal.stickyStruct('EVFlag',opts.EVFlag,...
        coder.internal.stickyStruct('Fs',opts.Fs,...
        coder.internal.stickyStruct('nfft',opts.nfft,...
        coder.internal.stickyStruct('p_eff',p_eff,...
        coder.internal.stickyStruct('eigenvals',eigenvals,...
        coder.internal.stickyStruct('signal_eigenvects',signal_eigenvects,...
        coder.internal.stickyStruct('noise_eigenvects',noise_eigenvects,...
        coder.internal.stickyStruct())))))))));
end
%-----------------------------------------------------------------------------------------
function [eigenvals,eigenvects] = computeeig(x,CorrFlag,CorrMatrOrd,nw,noverlap,window,EVFlag)
%COMPUTEEIG  Compute eigenvalues and eigenvectors of correlation matrix.
%
%   Inputs:
%      
%     x        - input vector or matrix
%     CorrFlag - (flag) indicates whether x is a correlation matrix
%     nw       - (integer) length of the rows of the data matrix 
%                (only used if x is vector)
%     noverlap - (integer) overlap between the rows of the data matrix
%                (used in conjunction with nw) 
%     window   - (vector) window to be applied to each column of data
%                 matrix  (not used if x is a correlation matrix)
%     EVFlag   -  True if eigenvector method, false if MUSIC.
%   
%
%   Outputs:
%    
%     eigenvals
%     eigenvects
%
%     If x is a matrix,  
%          If CorrFlag = 1, input x is a correlation matrix, we compute the
%          eigendecomposition and order the eigenvalues and eigenvectors.
%
%     If x is a vector,
%          a data matrix is formed by calling corrmtx unless a custom nw
%          and noverlap are specified. In that case, we use buffer to form
%          the data matrix. 
%
%     If window is not empty, each row of the data matrix will be
%     multiplied by the window.

coder.internal.prefer_const(CorrFlag,CorrMatrOrd,nw,noverlap,EVFlag);
isMATLAB = coder.target('MATLAB');
if CorrFlag && ~isvector(x)
    % Input is Correlation matrix
    coder.internal.assert(size(x,1) == size(x,2),...
                         'signal:music:CorrNeedsSquareMatrix')
    % Compute the eigenvectors and eigenvalues
    %[E,D] = eig((x+x')/2); % Ensure Hermitian
    % eig and schur are the same for hermitian inputs. In the generated 
    % C/C++ code, the outputs of eig are always complex. Using schur 
    % instead will produce real outputs for real inputs.
    
    [E,D]  = schur((x+x')/2);
    % D is always real in MATLAB as it contains the eigenvalues of a hermitan
    % matrix. However if input is complex, then the output in generated
    % C/C++ code will be complex with 0i as the imaginary part. Use real
    % cast here to ensure that eigenvals is real in the generated code.
    [eigenvals,indx] = sort(real(diag(D)),'descend');
    eigenvects = E(:,indx);

else
    if ~isvector(x)
        % Input is already a data matrix
        [Mx,Nx] = size(x); % Determine size of data matrix
        coder.internal.errorIf(EVFlag &&(Nx > Mx),...
            'signal:music:MatrixCannotHaveMoreColumsThanRows');
        y = x;
    else
        % x is a vector
        xcol = x(:); % Make it a column
        if coder.internal.isConst(isempty(nw)) && isempty(nw)
            y = corrmtx(xcol,CorrMatrOrd-1,'cov');
        else
            coder.internal.assert(length(xcol) > nw,'signal:music:invalidSegmentLength')
            coder.internal.errorIf(EVFlag && ...
                nw > (ceil((length(xcol)-nw)/(nw-noverlap))+1),'signal:music:invalidDataMatrix')
            Lx = length(xcol);
            y = buffNodelay(xcol,nw,noverlap)'./sqrt(Lx-nw); % Scale appropriately such that Y'*Y is a scaled estimate of R
        end
    end
    if coder.internal.isConst(isempty(window)) && isempty(window)
         dataMatrix = y;        
    else
        % Apply window to each row of data matrix
        coder.internal.assert(length(window) == size(y,2),...
                             'signal:music:InvalidDimensions');
        if isMATLAB                 
           dataMatrix = y.*repmat(window(:).',size(y,1),1);
        else
           eg = coder.internal.scalarEg(y,window);
           dataMatrix = coder.nullcopy(zeros(size(y),'like',eg));
           for j = 1:coder.internal.indexInt(size(y,2))
               for i = 1:coder.internal.indexInt(size(y,1))
                   dataMatrix(i,j) = y(i,j)*window(j);
               end
           end
        end       
    end

    % Compute the eigenvectors and eigenvalues via the SVD
    [~,S,eigenvects] = svd(dataMatrix,0);
    eigenvals = diag(S).^2; % We need to square the singular values here
end


%--------------------------------------------------------------------------------------------
function p_eff = determine_signal_space(p,eigenvals)
%DETERMINE_SIGNAL_SPACE   Determines the effective dimension of the signal subspace.
%   
%   Inputs:
%
%     p         - (scalar or vector) signal subspace dimension 
%                 (but may contain a desired threshold).
%     eigenvals - (vector) contains the eigenvalues (sorted in decreasing order)
%                 of the correlation matrix
%
%   Outputs:
%
%     p_eff - The effective dimension of the signal subspace. If a threshold
%             is given as p(2), the signal subspace will be equal to the number
%             of eigenvalues, NEIG, greater than the threshold times the smallest
%             eigenvalue. However, the dimension of the signal subspace is at most
%             p(1), so that if NEIG is greater than p(1), p_eff will be equal to
%             p(1). If the threshold criteria results in an empty signal subspace,
%             once again we make p_eff = p(1).


% Use the signal space dimension or the threshold to separate the noise subspace eigenvectors
if length(p) == 2
   % The threshold will be the input threshold times the smallest eigenvalue
   thresh = p(2)*eigenvals(end); 
   if coder.target('MATLAB')
       indx = find(eigenvals > thresh);
       if ~isempty(indx)
           p_eff = min( p(1), length(indx) );
       else
           p_eff = p(1);
       end
   else
       one = coder.internal.indexInt(1);
       indx = coder.internal.indexInt(0);
       for i = 1:coder.internal.indexInt(length(eigenvals))
           if eigenvals(i) > thresh
               indx = indx + one;
           end
       end
       if indx > 0
           p_eff = min(p(1),cast(indx,class(p)));
       else
           p_eff = p(1);
       end
   end
else
   p_eff = p(1);
end

function y = buffNodelay(x,nw,novp)
nwin = double(nw);
noverlap = double(novp);
if coder.target('MATLAB')
    y = buffer(x,nwin,noverlap,'nodelay');
else
    % equivalent implementation for y = buffer(x,nwin,noverlap,'nodelay');   
    coder.internal.prefer_const(nwin,noverlap);
    coder.internal.assert(nwin > noverlap,'signal:music:InvalidNoverlap')
    nx = length(x);
    h = nwin - noverlap;
    numCols = (nx-noverlap)/(nwin-noverlap);
    ytemp =signal.internal.stft.getSTFTColumns(x(:),nx,nwin,noverlap,1);% dummy value of 1 for Fs
    if numCols == floor(numCols)
        y = ytemp;
    else
        colEnd = floor(numCols);
        lastCol = zeros(nwin,1,'like',x);
        % copy the last noverlap number of samples from ytemp
        lastCol(1:noverlap) = ytemp(h+1:end,end);
        indexLast = nwin+h*(colEnd-1);%x(indexLast) is the last element of ytemp
        remainingElements = nx-indexLast;
        %copy remaining data in x.
        lastCol(noverlap+1:noverlap+remainingElements) = x(indexLast+1:end);
        y = [ytemp lastCol];
    end
end
% [EOF] - music.m
