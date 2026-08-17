function tf = keepOnlyVisualizationData(options)
%KEEPONLYVISUALIZATIONDATA True when only plot data should be retained.
%   The flag is intentionally opt-in so existing experiments keep their
%   diagnostic histories unless a script requests the lightweight mode.

tf = false;
if nargin < 1 || isempty(options) || ~isstruct(options)
    return;
end

if isTruthyScalar(options, 'keepOnlyVisualizationData') || ...
        isTruthyScalar(options, 'visualizationOnly')
    tf = true;
    return;
end

if isfield(options, 'historyMode') && isTextScalar(options.historyMode)
    mode = lower(strtrim(char(options.historyMode)));
    tf = any(strcmp(mode, {'visualization', 'visual', 'plot', 'minimal'}));
end
end

function tf = isTruthyScalar(s, fieldName)
tf = false;
if ~isfield(s, fieldName)
    return;
end

value = s.(fieldName);
tf = (islogical(value) || isnumeric(value)) && isscalar(value) && value ~= 0;
end

function tf = isTextScalar(value)
tf = ischar(value) || (isa(value, 'string') && isscalar(value));
end
