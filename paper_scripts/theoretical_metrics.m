function [MAE, RMSE, MCE] = theoretical_metrics(signal, prediction)
    errors = signal - prediction;
    MAE = mean(abs(errors));
    RMSE = sqrt(mean(errors.^2));
    MCE = mean(abs(errors).^3)^(1/3);
end
