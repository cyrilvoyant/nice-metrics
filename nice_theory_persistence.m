function p = nice_theory_persistence(sigma, rho)
    % nice_theory_persistence  Theoretical L^k errors of persistence for a weakly stationary series.
    %   p = nice_theory_persistence(sigma, rho)
    %   sigma : standard deviation of the series
    %   rho   : autocorrelation at the forecast lag (rho(1) for one step ahead)
    %
    % The persistence error D = y(t) - y(t - lag) has variance sigma_D^2 = 2 sigma^2 (1 - rho).
    % If D is Gaussian:
    %   MAE_P  = E|D|            = sqrt(2/pi) * sigma_D
    %   RMSE_P = sqrt(E[D^2])    = sigma_D
    %   RMCE_P = (E|D|^3)^(1/3)  = (2*sqrt(2)/sqrt(pi))^(1/3) * sigma_D
    % since E|D|^3 = 2*sqrt(2)/sqrt(pi) * sigma_D^3 for a centred Gaussian variable.
    %
    % Note: Eq. (10) and Eqs. (B.5)-(B.7) of the NICE^k article (author version arXiv:2508.01457v1) write
    % 4*sqrt(2)/sqrt(pi) for this moment; the factor is 2*sqrt(2)/sqrt(pi). See README.md.
    % The theoretical NICE^k of a model is then E_k(model) / p.E<k>.

    sigma_D = sqrt(2 * sigma .^ 2 .* (1 - rho));
    p.sigma_D = sigma_D;
    p.E1 = sqrt(2 / pi) * sigma_D;
    p.E2 = sigma_D;
    p.E3 = (2 * sqrt(2) / sqrt(pi)) ^ (1 / 3) * sigma_D;
    p.E3_article_constant = (4 * sqrt(2) / sqrt(pi)) ^ (1 / 3) * sigma_D; % as printed, for comparison only
end
