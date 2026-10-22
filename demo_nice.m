% demo_nice  Teaching demonstration of the NICE^k metrics on a synthetic solar-like signal.
%
% Part 1  One realisation: signal, forecasts of four models, classical scores and NICE^k.
% Part 2  Monte Carlo (200 realisations, random period and noise memory): theoretical versus
%         empirical errors of persistence, with the correct Gaussian constant for the cubic error
%         and with the constant printed in the article (see README.md, "Note on Eq. (B.5)").
%
% Outputs (folder results/): demo_series.csv, demo_scores.csv, montecarlo_persistence.csv.
% Figures: run  python make_figures.py  afterwards (or read the CSV files in any tool).
% Requires MATLAB R2021a or later; autocorr is replaced by a plain Pearson correlation.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'paper_scripts'));
out = fullfile(here, 'results');
if ~exist(out, 'dir'), mkdir(out); end

%% Part 1 - one realisation (same generator as the article)
rng(1);
T = 24; n_cycles = 500; filter_length = 3; h = 1;
[x, cs] = synthetic_signal(T, n_cycles, filter_length);
n = numel(x);

f.P      = predict_P(x, h);                           % simple persistence
f.SP     = predict_P_smart(x, cs, h);                 % smart persistence (additive clear-sky correction)
f.AR     = max(0, predict_AR(x, cs, h, 1));           % AR(1) on the clear-sky residual
f.Pcycle = predict_P_cyclic(x, h, T);                 % same hour of the previous cycle

valid = (T + 1:n)';                                   % common scored samples
names = fieldnames(f);
S = table();
for i = 1:numel(names)
    yh = f.(names{i});
    e = x(valid) - yh(valid);
    m = nice_metrics(x(valid), yh(valid), f.P(valid));
    S = [S; table(string(names{i}), ...
        sqrt(mean(e .^ 2)) / mean(x(valid)), mean(abs(e)) / mean(x(valid)), mean(e) / mean(x(valid)), ...
        1 - sum(e .^ 2) / sum((x(valid) - mean(x(valid))) .^ 2), ...
        m.NICE1, m.NICE2, m.NICE3, m.NICE_Sigma, ...
        'VariableNames', {'model', 'nRMSE', 'nMAE', 'nMBE', 'R2', 'NICE1', 'NICE2', 'NICE3', 'NICE_Sigma'})]; %#ok<AGROW>
end
disp(S);
writetable(S, fullfile(out, 'demo_scores.csv'));

days = (3 * T + 1:6 * T)';                            % three days for the figure
writetable(table(days, x(days), cs(days), f.P(days), f.SP(days), f.AR(days), f.Pcycle(days), ...
    'VariableNames', {'t', 'observed', 'clear_sky', 'P', 'SP', 'AR', 'Pcycle'}), fullfile(out, 'demo_series.csv'));

%% Part 2 - theory versus measurement for persistence, 200 realisations
rng(2);
n_runs = 200;
R = zeros(n_runs, 10);
for run = 1:n_runs
    Tr = randi([12, 96]);
    Lr = randi([2, max(2, round(Tr / 3))]);
    [xr, ~] = synthetic_signal(Tr, 300, Lr);
    s = xr(Tr + 1:end);                               % scored samples
    d = s(2:end) - s(1:end - 1);                      % persistence errors at lag 1
    rho1 = corr(s(1:end - 1), s(2:end));
    p = nice_theory_persistence(std(s), rho1);
    R(run, :) = [Tr, Lr, rho1, mean(abs(d)), sqrt(mean(d .^ 2)), mean(abs(d) .^ 3) ^ (1 / 3), ...
                 p.E1, p.E2, p.E3, p.E3_article_constant];
end
M = array2table(R, 'VariableNames', {'T', 'filter_length', 'rho1', 'MAE_P', 'RMSE_P', 'RMCE_P', ...
    'MAE_P_theory', 'RMSE_P_theory', 'RMCE_P_theory', 'RMCE_P_article_constant'});
writetable(M, fullfile(out, 'montecarlo_persistence.csv'));
fprintf('\nMedian ratio measured / theory over %d realisations:\n', n_runs);
fprintf('  MAE  %.3f\n  RMSE %.3f\n  RMCE %.3f (correct constant)   %.3f (constant as printed)\n', ...
    median(M.MAE_P ./ M.MAE_P_theory), median(M.RMSE_P ./ M.RMSE_P_theory), ...
    median(M.RMCE_P ./ M.RMCE_P_theory), median(M.RMCE_P ./ M.RMCE_P_article_constant));

function [x, cs] = synthetic_signal(T, n_cycles, filter_length)
    % Clear-sky-like cycle times bounded, smoothed Gaussian noise (generator of the article).
    t = (1:T * n_cycles)';
    w = 1 + randn(numel(t), 1);
    w = min(max(w, 0.2), 1.1);
    noise = filter(ones(filter_length, 1) / filter_length, 1, w);
    cs = max(0, 1000 * sin(2 * pi * t / T));
    x = max(0, cs .* noise);
end
