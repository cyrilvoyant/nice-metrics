function predictions = predict_AR(noisy_signal, signal, n_delta_t, p)
    % Modèle AR basé sur des boucles explicites.
    % Paramètres :
    % - noisy_signal : Série temporelle bruitée
    % - signal : Composante saisonnière
    % - n_delta_t : Horizon de prédiction
    % - p : Ordre du modèle AR
    %
    % Retour :
    % - predictions : Valeurs prédites avec recoloration

    % Taille du signal et calcul du résidu
    n_points = length(noisy_signal);
    residual_signal = noisy_signal - signal;
    average=mean(residual_signal);
    residual_signal=residual_signal-average;

    % Vérification de la condition sur p et n_delta_t
    if n_points - p - n_delta_t <= 0
        error('p + n_delta_t doit être inférieur à la longueur du signal.');
    end

    % Construction de X (matrice des retards) et y (valeurs cibles)
    X = zeros(n_points , p); % Matrice des retards (+1 pour la constante)
    y = zeros(n_points , 1);     % Vecteur des cibles

    for i = p:(n_points - n_delta_t)
        % Remplir chaque ligne de X avec les p valeurs précédentes
        for lag = 1:p
            X(i, lag) = residual_signal(i - lag+1);
        end

        % Remplir y avec la valeur à prédire (horizon n_delta_t)
        y(i) = residual_signal(i + n_delta_t);
    end

    % Estimation des coefficients AR par moindres carrés
   ar_coefficients=((X'*X)^-1)*X'*y;
 

    % Initialisation pour la prédiction des résidus
    predictions_residual = zeros(n_points, 1);
    for t =  1:(n_points - n_delta_t)
        % Calcul de la prédiction du résidu pour t + n_delta_t
        predictions_residual(t,1 ) = X(t,:)*ar_coefficients;
    end

    % Ajout de la composante saisonnière décalée
    predictions = zeros(n_points, 1);
    for t = 1:(n_points - n_delta_t)
        predictions(t + n_delta_t) = predictions_residual(t ) + average + signal(t + n_delta_t);
    end
end
