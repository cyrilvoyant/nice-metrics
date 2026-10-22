%Monte Carlo Simulation for NICE Metrics 
clear; clc; close all;
n_runs = 50;
metrics_empirical = struct();
predictor_labels = {'P', 'P_s', 'AR'};
metric_labels = {'NICE1', 'NICE2', 'NICE3', 'NICE_Sigma', 'nRMSE', 'nMAE', 'nMBE', 'R2'};
weights = [1/3, 1/3, 1/3]; % Poids pour NICE^Sigma

% Pré-allocation de la structure pour plus d'efficacité
for i = 1:length(predictor_labels)
    for j = 1:length(metric_labels)
        metrics_empirical.(predictor_labels{i}).(metric_labels{j}) = zeros(1, n_runs);
    end
end

for run = 1:n_runs
    % Génération des paramètres aléatoires
    T = randi([12, 240]); % Période aléatoire entre 12 et 240
    filter_length = randi([2, round(T/3)]); % Longueur du filtre entre T/4 et T/2
    
    % Génération du signal
    n_cycles = 500; % Nombre de cycles
    n_points = n_cycles * T;
    time = (1:n_points)';
    white_noise = 1 + randn(n_points, 1);
    white_noise(white_noise < 0.2) = 0.2;
    white_noise(white_noise > 1.1) = 1.1;
    correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
    signal = 1000 * sin(2 * pi * time / T); % Signal clair
    signal = max(0, signal); % Positivité garantie
    noisy_signal = signal .* correlated_noise; % Signal bruité
    noisy_signal = max(0, noisy_signal); % Positivité garantie
    
    % Prédictions des modèles
    predictions_P = predict_P(noisy_signal, 1); % Persistence simple
    predictions_P_smart = predict_P_smart(noisy_signal, signal, 1); % Smart persistence
    predictions_AR = predict_AR(noisy_signal, signal, 1, round(T/12)); % AR model
    predictions_AR = max(0, predictions_AR); 
    
    % Calcul des métriques de base
    start_index = max([T+1, 2]); % Indices valides pour ignorer les données initiales
    valid_noisy_signal = noisy_signal(start_index:end);
    mean_signal = mean(valid_noisy_signal);
    
    % Liste des prédictions
    predictors = {predictions_P, predictions_P_smart, predictions_AR};
    
    for pred_idx = 1:3
        predictions = predictors{pred_idx};
        
        % Ajustement des longueurs (CORRECTION IMPORTANTE)
        len_valid = min(length(valid_noisy_signal), length(predictions(start_index:end)));
        valid_signal_adjusted = valid_noisy_signal(1:len_valid);
        predictions_adjusted = predictions(start_index:start_index + len_valid - 1);
        
        % Calcul des erreurs
        errors = valid_signal_adjusted - predictions_adjusted; 
        
        % Calcul des métriques empiriques
        MAE_emp = mean(abs(errors));
        RMSE_emp = sqrt(mean(errors.^2));
        MCE_emp = mean(abs(errors).^3)^(1/3);
        MBE_emp = mean(errors); % Mean Bias Error
        nMAE_emp = MAE_emp / mean_signal; % Normalized MAE
        nRMSE_emp = RMSE_emp / mean_signal; % Normalized RMSE
        nMBE_emp = MBE_emp / mean_signal; % Normalized MBE
        
        % Calcul du coefficient de détermination R²
        ss_total = sum((valid_signal_adjusted - mean(valid_signal_adjusted)).^2);
        ss_residual = sum(errors.^2);
        if ss_total == 0 % Gérer le cas où la variance du signal est nulle
            R2_emp = NaN; % Ou une autre valeur appropriée, comme -Inf
        else
            R2_emp = 1 - (ss_residual / ss_total);
        end

        % Calcul des NICE empiriques
        if pred_idx == 1
            % Référence pour la persistance
            NICE1_emp = 1; NICE2_emp = 1; NICE3_emp = 1; NICE_Sigma_emp = 1;
        else
            if ~isfield(metrics_empirical.P, 'MAE') || isempty(metrics_empirical.P.MAE(run)) || metrics_empirical.P.MAE(run) == 0
                NICE1_emp = NaN;
            else
                NICE1_emp = MAE_emp / metrics_empirical.P.MAE(run);
            end
            if ~isfield(metrics_empirical.P, 'RMSE') || isempty(metrics_empirical.P.RMSE(run)) || metrics_empirical.P.RMSE(run) == 0
                NICE2_emp = NaN;
            else
                NICE2_emp = RMSE_emp / metrics_empirical.P.RMSE(run);
            end
            if ~isfield(metrics_empirical.P, 'MCE') || isempty(metrics_empirical.P.MCE(run)) || metrics_empirical.P.MCE(run) == 0
                NICE3_emp = NaN;
            else
                NICE3_emp = MCE_emp / metrics_empirical.P.MCE(run);
            end
            NICE_Sigma_emp = weights(1)*NICE1_emp + weights(2)*NICE2_emp + weights(3)*NICE3_emp;
        end
        
        % Stockage des métriques empiriques
        metrics_empirical.(predictor_labels{pred_idx}).MAE(run) = MAE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).RMSE(run) = RMSE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).MCE(run) = MCE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).MBE(run) = MBE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).nMAE(run) = nMAE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).nRMSE(run) = nRMSE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).nMBE(run) = nMBE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE1(run) = NICE1_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE2(run) = NICE2_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE3(run) = NICE3_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE_Sigma(run) = NICE_Sigma_emp;
        metrics_empirical.(predictor_labels{pred_idx}).R2(run) = R2_emp;
    end
end

% Tracé des prédictions AR vs Signal Bruité Valide
figure;
hold on;

% Vérifiez les tailles des séries
valid_signal = valid_noisy_signal; % Alias pour éviter la confusion
ar_predictions = predictions_AR(start_index:end);

% Subsample pour des tracés plus lisibles
subsample_idx = 1:1:min(length(valid_signal), length(ar_predictions));

% Tracé du signal bruité valide
plot(subsample_idx, valid_signal(subsample_idx), 'k', 'LineWidth', 1.5, 'DisplayName', 'Valid Noisy Signal');

% Tracé des prédictions AR
plot(subsample_idx, ar_predictions(subsample_idx), 'r--', 'LineWidth', 1.5, 'DisplayName', 'AR Prediction');

% Configuration des axes
xlabel('Time Index', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Signal Value', 'Interpreter', 'latex', 'FontSize', 14);
title('AR Predictions vs Valid Noisy Signal', 'Interpreter', 'latex', 'FontSize', 16);

% Légende
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'Best');

% Mise en forme
grid on;
set(gca, 'FontSize', 12, 'GridLineStyle', '--');
hold off;


%% Violin Plot with Kolmogorov-Smirnov Test for p-values
metrics_to_plot = {'NICE1', 'NICE2', 'NICE3', 'NICE_Sigma', 'nRMSE', 'nMAE', 'nMBE', 'R2'}; % Noms des champs dans la structure
latex_labels = {'NICE^1', 'NICE^2', 'NICE^3', 'NICE^{\Sigma}', 'nRMSE', 'nMAE', 'nMBE', 'R^2'}; % Labels en LaTeX
predictors_to_plot = {'P_s', 'AR'}; % Prédicteurs pour les violin plots
n_metrics = length(metrics_to_plot);
colors = [0.3 0.6 0.8; 0.9 0.3 0.4]; % Couleurs pour Ps (bleu) et AR (rouge)

% Création de la figure
figure;
hold on;

for i = 1:n_metrics
    metric = metrics_to_plot{i}; % Nom du champ

    % Extraire les données des deux prédicteurs
    data_P_s = metrics_empirical.(predictors_to_plot{1}).(metric)(:); % Données pour Ps
    data_AR = metrics_empirical.(predictors_to_plot{2}).(metric)(:); % Données pour AR

    % Calcul de la p-value avec le test Kolmogorov-Smirnov
    [~, p_value] = kstest2(data_P_s, data_AR); % Test KS pour comparer les distributions

    % Tracé des violins pour chaque prédicteur
    for j = 1:2
        predictor = predictors_to_plot{j};
        data = metrics_empirical.(predictor).(metric)(:); % Extraire les données
        pos = (i - 1) * 2 + j; % Position x pour les violins
        al_goodplot(data, pos, 0.5, colors(j, :), 'bilateral', [], [], 0.1);
    end

    % Ajout de la p-value au-dessus des violins (2 décimales)
    mid_point = (i - 1) * 2 + 1.5; % Point central entre les deux violins
    text(mid_point, 1.0, sprintf('$p=%.4f$', p_value), 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'center', 'FontSize', 18, 'Color', 'k');
end

% Configuration des axes
xticks(1.5:2:(n_metrics * 2 - 0.5));
xticklabels(latex_labels); % Utilisation des labels en LaTeX
xlabel('Error Metrics', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Values', 'Interpreter', 'latex', 'FontSize', 20);
title('Distribution of Normalized Error Metrics with $p$-values (KS Test)', 'Interpreter', 'latex', 'FontSize', 22);

% Légende propre
scatter(nan, nan, 40, 'MarkerFaceColor', colors(1, :), 'MarkerEdgeColor', 'k', 'DisplayName', '$P_s$');
scatter(nan, nan, 40, 'MarkerFaceColor', colors(2, :), 'MarkerEdgeColor', 'k', 'DisplayName', 'AR');
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'Best');

% Mise en forme
grid on;
set(gca, 'FontSize', 18, 'GridLineStyle', '--', 'XTickLabelRotation', 45);
xlim([0.5, n_metrics * 2 + 0.5]);
ylim([-0.1, 1.0]); % Ajustement de la hauteur pour inclure les p-values
hold off;
