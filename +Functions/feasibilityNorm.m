function value = feasibilityNorm(problem, x, y)
%FEASIBILITY_NORM Compute the split-constraint residual norm.
%   VALUE = FEASIBILITY_NORM(PROBLEM, X, Y) returns the Euclidean norm of the
%   current equality-constraint violation.
%
%   Input:
%       problem - problem structure with type 'LAD', 'lasso', 'elastic',
%                 or 'TV'.
%       x,y     - current split primal variables.
%
%   Output:
%       value   - ||Kx-y||_2 for LAD, ||Kx-y||_2 or ||x-y||_2 for LASSO,
%                 ||x-y||_2 for elastic-net, or
%                 ||D(x)-y||_2 for TV denoising.
	switch problem.type
	    case 'LAD'
	        value = norm(problem.K * x - y);
	    case 'lasso'
            if numel(y) == size(problem.K, 1)
	            value = norm(problem.K * x - y);
            elseif numel(y) == size(problem.K, 2)
	            value = norm(x - y);
            else
                error('feasibilityNorm:InvalidLassoSplit', ...
                    'LASSO y has incompatible dimension.');
            end
	    case 'elastic'
	        value = norm(x - y);
	    case 'TV'
	        value = norm(problem.D(x) - y);
	end
end
