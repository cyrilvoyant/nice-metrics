function m = nice_metrics(y, y_hat, y_ref, weights)
    % nice_metrics  Empirical NICE^k metrics of a forecast against a reference forecast.
    %   m = nice_metrics(y, y_hat)          reference = simple persistence y(t - 1) on the same samples
    %   m = nice_metrics(y, y_hat, y_ref)   any reference forecast aligned with y
    %   m = nice_metrics(..., weights)      weights of NICE^Sigma (default [1 1 1] / 3)
    %
    %   L^k error:  E_k = ( mean |y - y_hat|^k )^(1/k),  k = 1 (MAE), 2 (RMSE), 3 (RMCE)
    %   NICE^k    = E_k(model) / E_k(reference)
    %   NICE^Sigma = weighted sum of NICE^1, NICE^2 and NICE^3
    % Values below 1 beat the reference; 0 is a perfect forecast.
    % Samples where y, y_hat or y_ref is NaN are excluded from every metric.

    y = y(:); y_hat = y_hat(:);
    if nargin < 3 || isempty(y_ref)
        y_ref = [NaN; y(1:end - 1)];
    end
    y_ref = y_ref(:);
    if nargin < 4 || isempty(weights), weights = [1 1 1] / 3; end
    if numel(y_hat) ~= numel(y) || numel(y_ref) ~= numel(y)
        error('nice_metrics:size', 'y, y_hat and y_ref must have the same length.');
    end
    ok = ~isnan(y) & ~isnan(y_hat) & ~isnan(y_ref);
    e = y(ok) - y_hat(ok);
    r = y(ok) - y_ref(ok);
    for k = 1:3
        m.(sprintf('E%d_model', k)) = mean(abs(e) .^ k) ^ (1 / k);
        m.(sprintf('E%d_ref', k)) = mean(abs(r) .^ k) ^ (1 / k);
        m.(sprintf('NICE%d', k)) = m.(sprintf('E%d_model', k)) / m.(sprintf('E%d_ref', k));
    end
    m.NICE_Sigma = weights(1) * m.NICE1 + weights(2) * m.NICE2 + weights(3) * m.NICE3;
    m.n = nnz(ok);
end
