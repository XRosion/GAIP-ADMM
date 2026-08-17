function y = softThreshold(x, kappa)
%SOFT_THRESHOLD Apply the componentwise soft-thresholding operator.
%   Y = SOFT_THRESHOLD(X, KAPPA) computes
%       Y_i = sign(X_i) * max(abs(X_i) - KAPPA, 0)
%   for every component of X. This is the proximal operator of
%   KAPPA*||.||_1 and is used by the LAD and elastic-net subproblems.
%
%   Input:
%       x     - scalar, vector, or matrix of real values.
%       kappa - nonnegative threshold. It may be scalar or array-compatible
%               with x under MATLAB elementwise operations.
%
%   Output:
%       y     - soft-thresholded value with the same size as x.
%
%   Example:
%       softThreshold([-3; -0.5; 0; 2; 4], 1)
%       returns [-2; 0; 0; 1; 3].
y = sign(x) .* max(abs(x) - kappa, 0);
end
