function y = groupSoftThreshold(v, threshold)
%GROUPSOFTTHRESHOLD Block soft-thresholding for two-component TV gradients.
%   Y = GROUPSOFTTHRESHOLD(V, THRESHOLD) treats V as [v1; v2], where each
%   pair (v1_i, v2_i) is one isotropic-TV gradient block, and returns
%
%       max(1 - threshold / ||v_i||_2, 0) * v_i.

groupCount = numel(v) / 2;
if groupCount ~= floor(groupCount)
    error('groupSoftThreshold:InvalidLength', ...
        'Input length must be even for two-component groups.');
end

firstComponent = v(1:groupCount);
secondComponent = v(groupCount + 1:end);
groupMagnitude = sqrt(firstComponent.^2 + secondComponent.^2);
shrinkScale = max(1 - threshold ./ max(groupMagnitude, eps), 0);
y = [shrinkScale .* firstComponent; shrinkScale .* secondComponent];
end
