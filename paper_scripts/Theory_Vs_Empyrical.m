% MAIN SCRIPT - Monte Carlo Simulation for NICE Metrics
clear; clc; close all;

n_runs = 500;

metrics_theoretical = struct();
metrics_empirical = struct();

predictor_labels = {'P', 'SP', 'AR'};
metric_labels = {'NICE1', 'NICE2', 'NICE3', 'NICE_Sigma'};
weights = [1/3, 1/3, 1/3]; % Poids pour NICE^Sigma

for run = 1:n_runs
    % Génération des paramètres aléatoires
    T = randi([6, 240]); % Période aléatoire entre 6 et 240
    filter_length = randi([2, T]); % Longueur du filtre aléatoire entre 2 et T
    
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
    predictions_AR = predict_AR(noisy_signal,signal, 1, 1); % AR model
    
    % Calcul des métriques de base
    start_index = max([T+1, 2]); % Indices valides pour ignorer les données initiales
    rho = autocorr(noisy_signal(start_index:end), 'NumLags', 1);
    rho_1 = rho(2); % Autocorrélation au lag 1
    sigma_signal = std(noisy_signal(start_index:end)); % Écart-type du signal

    % Calcul des métriques théoriques pour la persistance
    MAE_P_th = sqrt(2/pi) * sqrt(2 * sigma_signal^2 * (1 - rho_1));
    RMSE_P_th = sqrt(2 * sigma_signal^2 * (1 - rho_1));
    MCE_P_th = (4*sqrt(2)/sqrt(pi) * (2 * sigma_signal^2 * (1 - rho_1))^(3/2))^(1/3);

    % Stockage des métriques théoriques pour la persistance
    metrics_theoretical.P.MAE(run) = MAE_P_th;
    metrics_theoretical.P.RMSE(run) = RMSE_P_th;
    metrics_theoretical.P.MCE(run) = MCE_P_th;

    % Liste des prédictions
    predictors = {predictions_P, predictions_P_smart, predictions_AR};
    
    for pred_idx = 1:3
        predictions = predictors{pred_idx};
        errors = noisy_signal(start_index:end) - predictions(start_index:end); % Erreurs
        
        % Calcul des métriques empiriques
        MAE_emp = mean(abs(errors));
        RMSE_emp = sqrt(mean(errors.^2));
        MCE_emp = mean(abs(errors).^3)^(1/3);
        
        % Calcul du coefficient R²
        ss_total = sum((noisy_signal(start_index:end) - mean(noisy_signal(start_index:end))).^2);
        ss_residual = sum(errors.^2);
        R2 = 1 - (ss_residual / ss_total);
        
        % Calcul des métriques théoriques pour les prédicteurs (si différent de P)
        if pred_idx ~= 1
            [MAE_pred_th, RMSE_pred_th, MCE_pred_th] = theoretical_metrics(noisy_signal(start_index:end), predictions(start_index:end));
            metrics_theoretical.(predictor_labels{pred_idx}).MAE(run) = MAE_pred_th;
            metrics_theoretical.(predictor_labels{pred_idx}).RMSE(run) = RMSE_pred_th;
            metrics_theoretical.(predictor_labels{pred_idx}).MCE(run) = MCE_pred_th;
        else
            % La persistance utilise ses propres métriques théoriques
            metrics_theoretical.(predictor_labels{pred_idx}).MAE(run) = MAE_P_th;
            metrics_theoretical.(predictor_labels{pred_idx}).RMSE(run) = RMSE_P_th;
            metrics_theoretical.(predictor_labels{pred_idx}).MCE(run) = MCE_P_th;
        end

        % Calcul des NICE (théorique)
        NICE1_th = metrics_theoretical.(predictor_labels{pred_idx}).MAE(run) / MAE_P_th;
        NICE2_th = metrics_theoretical.(predictor_labels{pred_idx}).RMSE(run) / RMSE_P_th;
        NICE3_th = metrics_theoretical.(predictor_labels{pred_idx}).MCE(run) / MCE_P_th;
        NICE_Sigma_th = weights(1)*NICE1_th + weights(2)*NICE2_th + weights(3)*NICE3_th;
        
        % Calcul des NICE (empirique)
        if pred_idx == 1
            % Référence pour la persistance
            NICE1_emp = 1; NICE2_emp = 1; NICE3_emp = 1; NICE_Sigma_emp = 1;
        else
            NICE1_emp = MAE_emp / metrics_empirical.P.MAE(run);
            NICE2_emp = RMSE_emp / metrics_empirical.P.RMSE(run);
            NICE3_emp = MCE_emp / metrics_empirical.P.MCE(run);
            NICE_Sigma_emp = weights(1)*NICE1_emp + weights(2)*NICE2_emp + weights(3)*NICE3_emp;
        end
        
        % Stockage des métriques empiriques
        metrics_empirical.(predictor_labels{pred_idx}).MAE(run) = MAE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).RMSE(run) = RMSE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).MCE(run) = MCE_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE1(run) = NICE1_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE2(run) = NICE2_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE3(run) = NICE3_emp;
        metrics_empirical.(predictor_labels{pred_idx}).NICE_Sigma(run) = NICE_Sigma_emp;
        metrics_empirical.(predictor_labels{pred_idx}).R2(run) = R2;

        % Stockage des métriques théoriques
        metrics_theoretical.(predictor_labels{pred_idx}).NICE1(run) = NICE1_th;
        metrics_theoretical.(predictor_labels{pred_idx}).NICE2(run) = NICE2_th;
        metrics_theoretical.(predictor_labels{pred_idx}).NICE3(run) = NICE3_th;
        metrics_theoretical.(predictor_labels{pred_idx}).NICE_Sigma(run) = NICE_Sigma_th;
    end
end

% Visualisation Améliorée avec R² colorés
figure;
titles = {'$\mathtt{NICE^1}$ ', ...
          '$\mathtt{NICE^2}$ ', ...
          '$\mathtt{NICE^3}$ ', ...
          '$\mathtt{NICE^{\Sigma}}$ '};

% Couleurs pour chaque prédicteur
colors = {[0, 0, 1, 0.5], [1, 0, 0, 0.5]}; % Bleu pour P_smart, Rouge pour AR

for i = 1:4
    subplot(2, 2, i);
    hold on;

    for pred_idx = 2:3 % P_smart (bleu) et AR (rouge)
        x = metrics_theoretical.(predictor_labels{pred_idx}).(metric_labels{i});
        y = metrics_empirical.(predictor_labels{pred_idx}).(metric_labels{i});
        
        % Calcul du R²
        R2 = 1 - sum((y - x).^2) / sum((y - mean(y)).^2);

        % Scatter plot
        scatter(x, y, 10, 'filled', 'MarkerFaceColor', colors{pred_idx - 1}(1:3), ...
                'MarkerFaceAlpha',0.5,'DisplayName', ['$\mathrm{' predictor_labels{pred_idx} '}$']);
        
        % Texte pour R²
        text(0.22, 0.92 - 0.2 * (pred_idx - 2), ...
             ['$\mathrm{R^2=' sprintf('%.3f', R2) '}$'], ...
             'Units', 'normalized', 'FontSize', 16, ...
             'Interpreter', 'latex', 'Color', colors{pred_idx - 1}(1:3));
    end

    % Droite y=x pour référence
    plot([0.2, 1], [0.2, 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'y=x');
    
    % Titres et axes en LaTeX
    title(titles{i}, 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Theoretical', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Empirical', 'Interpreter', 'latex', 'FontSize', 16);

    % Ajuster les limites des axes
    xlim([0.2, 1]);
    ylim([0.2, 1]);
    
    % Style de la grille
    grid on;
    set(gca, 'GridLineStyle', ':', 'FontSize', 14, 'FontName', 'Times');
    axis square;
end

% Légende ajustée
legend('Location', 'Best', 'Interpreter', 'latex', 'FontSize', 14);
% Sauvegarder la figure en PNG
saveas(gcf, 'NICE_Metrics_Comparison.png');
