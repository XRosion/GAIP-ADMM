function value = tvGroupNormSum(y)
%TVGROUPNORMSUM Sum isotropic-TV two-component block norms.
%   VALUE = TVGROUPNORMSUM(Y) evaluates sum_i sqrt(y1_i^2 + y2_i^2)
%   for a stacked vector Y = [y1; y2].

groupCount = numel(y) / 2;
if groupCount ~= floor(groupCount)
    error('tvGroupNormSum:InvalidLength', ...
        'Input length must be even for two-component TV groups.');
end

firstComponent = y(1:groupCount);
secondComponent = y(groupCount + 1:end);
value = sum(sqrt(firstComponent.^2 + secondComponent.^2));
end
