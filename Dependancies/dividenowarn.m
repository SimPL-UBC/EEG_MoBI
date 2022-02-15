function y = dividenowarn(num,den)
% DIVIDENOWARN Divides two polynomials while suppressing warnings.
% DIVIDENOWARN(NUM,DEN) array divides two polynomials but suppresses warnings 
% to avoid "Divide by zero" warnings.

%   Copyright 1988-2020 The MathWorks, Inc.
%#codegen
isMATLAB = coder.target('MATLAB');
if isMATLAB
   s = warning; % Cache warning state
   warning off  % Avoid "Divide by zero" warnings
   y = (num./den);
   warning(s);  % Reset warning state
else
   y = (num./den);
end
% [EOF] dividenowarn.m
