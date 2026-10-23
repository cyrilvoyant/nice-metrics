# NICE metrics: scientific use contract

## Purpose and authority

Evaluate a deterministic forecast against an explicitly identified reference on aligned samples. Read `nice_metrics.m`, `README.md` and `CITATION.cff` first. This document describes the empirical function at source commit `913c880a7d582d4a4365dd62a9b97e78f2ad359f`; it does not change the code. `websemantic.json` is a documentation descriptor, not a deployed WebSemantic adapter or an executable validation gate.

## Inputs and units

- `y`: observations, real numeric vector.
- `y_hat`: model predictions, same length, timestamps, forecast horizon and physical unit as `y`.
- `y_ref`: reference predictions aligned with those same observations. Ask which reference and horizon are intended. If omitted, the MATLAB code constructs `[NaN; y(1:end-1)]`: one sample lag, not automatically one hour or 24 hours. The first sample is then excluded.
- `weights`: three coefficients for `NICE_Sigma`. The code default is `[1 1 1]/3`; it neither normalizes nor validates supplied weights. For a convex summary, explicitly validate finite nonnegative weights summing to one before calling it.

The function is unit-agnostic. For solar irradiance use the same unit, for example W/m², in all three vectors. Do not silently mix power, energy, irradiance, normalized signals or time resolutions. Record the unit and forecast horizon with the result. Apply any day/night mask identically to observations and both forecasts.

## Equations and outputs

For k = 1, 2, 3, `E_k = mean(abs(y - forecast).^k)^(1/k)`. The output structure contains `E1_model` through `E3_model`, `E1_ref` through `E3_ref`, `NICE1` through `NICE3`, `NICE_Sigma`, and `n`.

`NICEk = Ek_model / Ek_ref`; `NICE_Sigma = sum(weights .* [NICE1 NICE2 NICE3])`. Error norms retain the input unit; NICE values and the weights are dimensionless; `n` is the count of retained samples. Below one means a smaller empirical error than the chosen reference for that norm. It is not a universal deployment recommendation.

## Missing data and limits

The implementation excludes a sample when any of `y`, `y_hat`, `y_ref` is NaN. It does not reject infinity, complex values or an empty retained set. Validate those cases in the caller. A zero reference error makes the ratio undefined or infinite; report that explicitly, never replace the denominator silently. The implementation does not estimate confidence intervals, prove significance, select a forecast horizon or prevent evaluation leakage. Compare models on the same retained samples.

The empirical function uses measured reference errors. `nice_theory_persistence.m` is a separate theoretical calculation with additional assumptions. See the README note about the Gaussian third absolute moment in the author manuscript: do not transfer a theoretical constant issue to the empirical ratios, and do not claim the publisher PDF has been checked if only the author version was read. Original `paper_scripts/` remain unchanged.

## Minimal example

From the repository directory, in MATLAB:

```matlab
m = nice_metrics([500 450 600 550], [460 480 560 580], [400 500 450 600]);
disp([m.NICE1 m.NICE2 m.NICE3 m.NICE_Sigma m.n]);
```

These are synthetic teaching irradiances in W/m², not field validation data. Expected values, rounded: 0.400, 0.365, 0.340, 0.368, 4. An agent must report whether it actually executed MATLAB; reading these values is not execution.

## Attribution

Article: C. Voyant, M. Despotovic, L. Garcia-Gutierrez, R. Amaro e Silva, P. Lauret, T. Soubdhan and N. Bailek, *NICEᵏ metrics: Unified and multidimensional framework for evaluating deterministic solar forecasting accuracy*, Sustainable Energy Technologies and Assessments 83 (2025), 104588. https://doi.org/10.1016/j.seta.2025.104588

Software concept DOI: https://doi.org/10.5281/zenodo.23562549 ; archived v1.0.0: https://doi.org/10.5281/zenodo.23562550 . Documentation added after that release is not automatically part of its archive. Cite the exact software commit or release used. License: MIT.
