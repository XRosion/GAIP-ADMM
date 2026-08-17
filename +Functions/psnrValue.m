function value = psnrValue(problem, x)
%PSNRVALUE Compute PSNR for a recovered image vector.
%   VALUE = PSNRVALUE(PROBLEM, X) uses PROBLEM.xtrue as the reference image.
%   Images are assumed to be scaled to [0, 1] unless PROBLEM.pixelMax is set.

if ~isfield(problem, 'xtrue')
    error('psnrValue:MissingReference', 'problem.xtrue is required.');
end

if isfield(problem, 'pixelMax') && ~isempty(problem.pixelMax)
    pixelMax = problem.pixelMax;
else
    pixelMax = 1;
end

mse = mean((x - problem.xtrue).^2);
value = 10 * log10(pixelMax^2 / max(mse, eps));
end
