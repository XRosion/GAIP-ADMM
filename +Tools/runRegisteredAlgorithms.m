function [results, elapsedTimes] = runRegisteredAlgorithms(Algorithms, problem, options, statusPrinter)
%RUNREGISTEREDALGORITHMS Run algorithms listed in a simple registry.
%   The registry follows the same convention used by Main*.m:
%       Type = 0: comparison algorithm, out = Func(problem, options)
%       Type = 1: proposed algorithm,   out = Func(problem, options, Parameters)

if nargin < 4
    statusPrinter = [];
end

numAlgorithms = numel(Algorithms);
results = cell(numAlgorithms, 1);
elapsedTimes = zeros(numAlgorithms, 1);

for iAlg = 1:numAlgorithms
    alg = Algorithms(iAlg);
    fprintf('\n>>> Algorithm %d/%d : %s (%s)\n', ...
        iAlg, numAlgorithms, alg.Name, alg.Description);

    tic;
    if alg.Type == 1
        if ~isfield(alg, 'Parameters') || isempty(alg.Parameters)
            error('runRegisteredAlgorithms:MissingParameters', ...
                'Algorithm "%s" has Type = 1 but no Parameters field.', alg.Name);
        end
        out = alg.Func(problem, options, alg.Parameters);
    else
        out = alg.Func(problem, options);
    end
    elapsedTimes(iAlg) = toc;

    out.name = alg.Name;
    out.displayName = sprintf('%s (%s)', alg.Name, alg.Description);
    if isfield(out, 'hist') && ~isfield(out, 'history')
        out.history = out.hist;
    elseif isfield(out, 'history') && ~isfield(out, 'hist')
        out.hist = out.history;
    end

    if ~isempty(statusPrinter)
        statusPrinter(out, elapsedTimes(iAlg), alg);
    end

    results{iAlg} = out;
end
end
