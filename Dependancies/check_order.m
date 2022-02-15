function [n_out, w, trivalwin] = check_order(n_in)
%CHECK_ORDER Checks the order passed to the window functions.
% [N,W,TRIVALWIN] = CHECK_ORDER(N_ESTIMATE) will round N_ESTIMATE to the
% nearest integer if it is not already an integer. In special cases (N is
% [], 0, or 1), TRIVALWIN will be set to flag that W has been modified.

%   Copyright 1988-2020 The MathWorks, Inc.

%#codegen

w = 0;
trivalwin = 0;

% Special case of N is []
if isempty(n_in)
    n_out = zeros(class(n_in));
    w = zeros(0,1);
    trivalwin = 1;
    return
end

validateattributes(n_in,{'numeric'},{'scalar','finite','real','nonnegative'},'check_order','N');
n_in = n_in(1);

% Check if order is already an integer or empty
% If not, round to nearest integer.
if n_in == floor(n_in)
    n_out = n_in;
else
    n_out = round(n_in);
    coder.internal.warning('signal:check_order:InvalidOrderRounding');
end
    
% special cases: N is 0 0r 1
if n_out == 0 
    w = zeros(0,1);       % Empty matrix: 0-by-1
    trivalwin = 1;
elseif n_out == 1
    w = 1;
    trivalwin = 1;
end


% LocalWords:  TRIVALWIN
