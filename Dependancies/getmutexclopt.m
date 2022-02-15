function [opt, arglist] = getmutexclopt(validopts,defaultopt,arglist)
%GETMUTEXCLOPT - get any of the specified mutually exclusive options and
% remove from the argument list.  Allows initial matches.
%
%   This function is for internal purposes only and may be removed in a
%   future release.
%
%   validtypes  - a cell array of valid options
%                 (e.g. {'power','ms','psd'})
%
%   defaulttype - the default option to use if no type is found
%
%   arglist    - the input argument list
%
%   Errors out if different estimation types are matched in the arglist.    
%
%   See also CHKUNUSEDOPT.

%   Copyright 2015 The MathWorks, Inc.

opt = defaultopt;
found = false;

iarg = 1;
while iarg <= numel(arglist)
  arg = arglist{iarg};
  if ischar(arg) && isrow(arg)
    matches = find(strncmpi(arg,validopts,length(arg)));
    if ~isempty(matches)
      if ~found
        found = true;
        opt = validopts{matches(1)};
        arglist(iarg) = [];
      else
        error(message('signal:getmutexclopt:ConflictingOptions', ...
                      opt,validopts{matches(1)}));
      end
    else
      iarg = iarg + 1;
    end
  else
   iarg = iarg + 1;
  end
end
    
    
    