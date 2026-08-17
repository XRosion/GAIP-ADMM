function tf = keepOnlyFinalVisualizationData(options)
%KEEPONLYFINALVISUALIZATIONDATA True when only final plot values are needed.

tf = Tools.keepOnlyVisualizationData(options) && ...
    isTruthyScalar(options, 'keepOnlyFinalVisualizationData');
end

function tf = isTruthyScalar(s, fieldName)
tf = false;
if nargin < 1 || isempty(s) || ~isstruct(s) || ~isfield(s, fieldName)
    return;
end

value = s.(fieldName);
if islogical(value) && isscalar(value)
    tf = value;
elseif isnumeric(value) && isscalar(value) && isfinite(value)
    tf = value ~= 0;
end
end
