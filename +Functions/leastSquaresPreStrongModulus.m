function modulus = leastSquaresPreStrongModulus(K)
%LEASTSQUARESPRESTRONGMODULUS Conservative modulus of 0.5*||K*x-b||^2.
%   The exact modulus is lambda_min(K'*K) = sigma_min(K)^2. A wide or
%   numerically rank-deficient matrix has modulus zero. For a full-column-
%   rank matrix, subtracting the rank tolerance from sigma_min gives a
%   conservative value for use in penalty-parameter updates.

validateattributes(K, {'numeric'}, {'2d', 'nonempty', 'real', 'finite'}, ...
    mfilename, 'K');

[rowCount, columnCount] = size(K);
if rowCount < columnCount
    modulus = 0;
    return;
end

singularValues = svd(K, 'econ');
largestSingularValue = singularValues(1);
smallestSingularValue = singularValues(end);
rankTolerance = max(size(K)) * eps(max(largestSingularValue, 1));
safeSmallestSingularValue = max(smallestSingularValue - rankTolerance, 0);
modulus = safeSmallestSingularValue^2;
end
