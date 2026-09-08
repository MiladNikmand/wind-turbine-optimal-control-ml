function out = available(models)
%PREDICTORS.AVAILABLE  Which predictors can this MATLAB install actually run?
%
%   predictors.available()                    print a table for all models
%   T = predictors.available()                return it
%   ok = predictors.available({'rbf','gp'})   logical vector
%
% Checks the toolboxes each predictor needs. Worth running once before a
% long batch: the alternative is discovering three hours in that ARIMA
% needs Econometrics, or that bilstm calls trainNetwork, which Mathworks
% REMOVED in R2024a.

    all_names = fieldnames(predictors.registry());
    if nargin == 0 || isempty(models)
        models = all_names;
        want_table = true;
    else
        if ischar(models), models = {models}; end
        want_table = false;
    end

    n = numel(models);
    name = cell(n,1); req = cell(n,1); ok = false(n,1); why = cell(n,1);

    for i = 1:n
        e = predictors.registry(models{i});
        name{i} = e.name;
        req{i}  = strjoin(e.requires, '+');
        if isempty(req{i}), req{i} = '-'; end

        [ok(i), why{i}] = local_check(e.requires);
    end

    if ~want_table && nargout > 0
        out = ok;
        return;
    end

    out = table(name, req, ok, why, ...
        'VariableNames', {'Model','Requires','Available','Reason'});

    if nargout == 0
        fprintf('\nPredictor availability on this MATLAB install:\n\n');
        disp(out);
        n_bad = sum(~ok);
        if n_bad > 0
            fprintf(2, '%d model(s) unavailable -- exclude them from models_to_run.\n\n', n_bad);
        else
            fprintf('All %d models available.\n\n', n);
        end
        clear out;
    end
end

% =========================================================================
function [ok, why] = local_check(reqs)
    ok = true; why = 'ok';
    for r = reqs
        switch r{1}
            case 'stats'
                if isempty(ver('stats'))
                    ok = false; why = 'needs Statistics and Machine Learning Toolbox'; return;
                end
            case 'econ'
                if isempty(ver('econ'))
                    ok = false; why = 'needs Econometrics Toolbox'; return;
                end
            case 'deep_legacy'
                if isempty(ver('nnet'))
                    ok = false; why = 'needs Deep Learning Toolbox'; return;
                end
                if exist('trainNetwork', 'file') ~= 2
                    ok = false; why = 'trainNetwork removed in R2024a+ -- use gru/lstm/tcn instead'; return;
                end
        end
    end
end
