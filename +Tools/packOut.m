function out = packOut(name, x, y, lambda, history, options)
%PACK_OUTPUT Standardize algorithm result structures.
%   OUT = PACK_OUTPUT(NAME, X, Y, LAMBDA, HISTORY) packages final iterates
%   and the recorded history in the common format consumed by plotting tools.
%
%   Input:
%       name   - algorithm display name.
%       x,y    - final primal split variables.
%       lambda - final multiplier estimate.
%       history - history structure assembled by appendHistory.
%
%   Output:
%       out    - structure with fields name, x, y, lambda, history, and iter.
if nargin < 6
    options = struct();
end

out.name = name;
if Tools.keepOnlyVisualizationData(options)
    out.history = history;
    out.iter = numel(history.obj);
    return;
end

out.x = x;
out.y = y;
out.lambda = lambda;
out.history = history;
out.iter = numel(history.obj);
end
