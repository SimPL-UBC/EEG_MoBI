function chknostropts(varargin)
%CHKNOSTROPTS Checks that no strings are specified
%   errors on the first encountered string
%
%   This file is for internal use only

%   Copyright 2013-2019 The MathWorks, Inc.

%#codegen

for i = 1:numel(varargin)
    coder.internal.errorIf(ischar(varargin{i}),...
        'signal:chknostropts:UnknownOption',varargin{i});
end