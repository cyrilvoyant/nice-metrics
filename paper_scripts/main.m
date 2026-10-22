% MAIN SCRIPT
clear;
clc;
close all;

% Paramètres
T = 24;         % Période (en heures)
n_cycles = 500; % Nombre de cycles
n_points = n_cycles * T; % Nombre total de points
n_delta_t = 1;  % Décalage temporel pour la corrélation cyclique
filter_length = 3; % Longueur du filtre de lissage
p = 1; % Ordre du modèle PAR(p) (1 * T)

% Générer le signal
time = (1:n_points)';

% Génération du bruit borné entre 0.2 et 1.1
white_noise = 1 + randn(n_points, 1); % Bruit gaussien centré autour de 1
white_noise(white_noise < 0.2) = 0.2; % Limite inférieure
white_noise(white_noise > 1.1) = 1.1; % Limite supérieure
correlated_noise = filter(ones(filter_length, 1) / filter_length, 1, white_noise);
signal = 1000 * sin(2 * pi * time / T); % Composante sinusoïdale (clear-sky)
signal = max(0, signal); % Assurer que signal est positif ou nul
noisy_signal = signal .* correlated_noise; % Signal final bruité
noisy_signal = max(0, noisy_signal); % Assurer que noisy_signal est positif ou nul

% Prédictions avec différents modèles
predictions_P = predict_P(noisy_signal, n_delta_t); % Persistence simple
predictions_P_smart = predict_P_smart(noisy_signal, signal, n_delta_t); % Smart persistence
predictions_PAR = predict_PAR(noisy_signal, T, n_delta_t, p); % Modèle PAR
predictions_P_cyclic = predict_P_cyclic(noisy_signal, n_delta_t, T); % Modèle cyclique

% Gestion des valeurs négatives pour les prédictions
predictions_P = max(0, predictions_P);
predictions_P_smart = max(0, predictions_P_smart);
predictions_PAR = max(0, predictions_PAR);
predictions_P_cyclic = max(0, predictions_P_cyclic);

% Calculer les métriques d'évaluation
start_index = max([T+1, n_delta_t+1]); % Ignorer les indices initiaux non valides
errors = zeros(length(noisy_signal), 4); % Matrice pour stocker les erreurs
metrics = struct();

% Calculer les erreurs pour chaque modèle
errors(:, 1) = noisy_signal - predictions_P;
errors(:, 2) = noisy_signal - predictions_P_smart;
errors(:, 3) = noisy_signal - predictions_PAR;
errors(:, 4) = noisy_signal - predictions_P_cyclic;

% Calculer les différences du signal bruité
diff_noisy_signal = diff(noisy_signal);

% Calculer rho (autocorrélation au lag 1) sur les différences (EMPIRIQUEMENT)
rho = autocorr(noisy_signal(start_index-1:end), 'NumLags', 1); % Correction ici
rho_1 = rho(2);

% Écart-type du SIGNAL (calculé EMPIRIQUEMENT)
sigma_signal = std(noisy_signal(start_index:end));

% Constantes de normalisation pour NICE théorique (CORRECTES selon la théorie fournie)
MAE_P_th = sqrt(2/pi) * sqrt(2 * sigma_signal^2 * (1 - rho_1));
RMSE_P_th = sqrt(2 * sigma_signal^2 * (1 - rho_1));
MCE_P_th = (4*sqrt(2)/sqrt(pi) * (2 * sigma_signal^2 * (1 - rho_1))^(3/2))^(1/3);

% Calcul des métriques
weights = [1/3, 1/3, 1/3]; % Poids pour NICE^Sigma
for i = 1:4
    % RMSE et nRMSE
    metrics.RMSE(i) = sqrt(mean(errors(start_index:end, i).^2));
    metrics.nRMSE(i) = metrics.RMSE(i) / mean(noisy_signal(start_index:end));
    
    % MAE et nMAE
    metrics.MAE(i) = mean(abs(errors(start_index:end, i)));
    metrics.nMAE(i) = metrics.MAE(i) / mean(noisy_signal(start_index:end));
    
    % MCE (Mean Cubic Error) - Calcul corrigé avec racine cubique
    metrics.MCE(i) = mean(abs(errors(start_index:end, i)).^3)^(1/3);
    
    % MBE et nMBE
    metrics.MBE(i) = mean(errors(start_index:end, i));
    metrics.nMBE(i) = metrics.MBE(i) / mean(noisy_signal(start_index:end));
    
    % R² (Coefficient de détermination)
    ss_total = sum((noisy_signal(start_index:end) - mean(noisy_signal(start_index:end))).^2);
    ss_residual = sum(errors(start_index:end, i).^2);
    metrics.R2(i) = 1 - (ss_residual / ss_total);

    % Calcul des métriques NICE théoriques (CORRECTES selon la théorie fournie)
    metrics.NICE1_theory(i) = metrics.MAE(i) / MAE_P_th;
    metrics.NICE2_theory(i) = metrics.RMSE(i) / RMSE_P_th;
    metrics.NICE3_theory(i) = metrics.MCE(i) / MCE_P_th;
end

% Calculs empiriques pour NICE
RMSE_P = metrics.RMSE(1); % Référence : Persistence simple
for i = 1:4
    % NICE^2 empirique
    metrics.NICE2_empirical(i) = metrics.RMSE(i) / RMSE_P; 
    
    % NICE^1 empirique
    metrics.NICE1_empirical(i) = metrics.MAE(i) / metrics.MAE(1);
    
    % NICE^3 empirique
    metrics.NICE3_empirical(i) = metrics.MCE(i) / metrics.MCE(1);
end

% NICE^Sigma: Empirique
metrics.NICE_Sigma_empirical = zeros(1, 4);
for i = 1:4
    metrics.NICE_Sigma_empirical(i) = weights(1) * metrics.NICE1_empirical(i) + ...
                                      weights(2) * metrics.NICE2_empirical(i) + ...
                                      weights(3) * metrics.NICE3_empirical(i);
end

% NICE^Sigma: Théorique
metrics.NICE_Sigma_theory = zeros(1, 4);
for i = 1:4
    metrics.NICE_Sigma_theory(i) = weights(1) * metrics.NICE1_theory(i) + ...
                                   weights(2) * metrics.NICE2_theory(i) + ...
                                   weights(3) * metrics.NICE3_theory(i);
end

% Affichage des résultats
disp('--- Metrics Summary ---');
disp(metrics);

% Visualisation des signaux et des prédictions
figure;
plot(time, noisy_signal, 'k', 'DisplayName', 'Noisy Signal');
hold on;
plot(time, predictions_P, 'b', 'DisplayName', 'Persistence');
plot(time, predictions_P_smart, 'r', 'DisplayName', 'Smart Persistence');
plot(time, predictions_PAR, 'g', 'DisplayName', 'PAR Model');
plot(time, predictions_P_cyclic, 'm', 'DisplayName', 'Cyclic Model');
xlabel('Time');
ylabel('Signal Value');
title('Forecasting Models and Noisy Signal');
legend;
grid on;
