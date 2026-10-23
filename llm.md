# NICE metrics: scientific use contract

## Purpose and authority

Evaluate a deterministic forecast against **simple persistence at the same forecast horizon**. This reference is part of the NICE definition, not a selectable benchmark. Read the article first (arXiv:2508.01457, Sections 2.1–2.5, especially Eqs. 2–6 and the horizon extension), then `README.md`, `CITATION.cff` and the implementation. This document describes the empirical function at source commit `913c880a7d582d4a4365dd62a9b97e78f2ad359f`; it does not change the code. `websemantic.json` is a documentation descriptor, not a deployed WebSemantic adapter or an executable validation gate.

The scientific motivation is to evaluate whether modeling adds value beyond repeating the last available observation, using a fixed data-driven reference and consistent Lk norms. A ratio against smart persistence, climatology, a cyclic forecast or any other chosen benchmark must not be presented as NICE as defined in this article. The generic `y_ref` argument permits such calculations mechanically; that flexibility does not redefine the method. Do not relabel `1-NICEk` as the article's metric or as a universal skill score.

## Inputs and units

- `y`: observations, real numeric vector.
- `y_hat`: model predictions, same length, timestamps, forecast horizon and physical unit as `y`.
- `y_ref`: **simple-persistence predictions only**, aligned with those observations. At horizon h samples, the prediction for target t is the observation at t-h. Construct it from the full regularly sampled series before applying the common evaluation mask. Ask for the horizon and sampling interval, not for a choice of benchmark. If omitted, the MATLAB code constructs `[NaN; y(1:end-1)]`: valid for one-step persistence only when the supplied y retains the original consecutive time grid. It is not automatically one hour or 24 hours. The first sample is then excluded. For h>1 or a preselected evaluation subset, supply the correctly aligned persistence vector explicitly.
- `weights`: three coefficients for `NICE_Sigma`. The code default is `[1 1 1]/3`; it neither normalizes nor validates supplied weights. The article requires nonnegative weights summing to one (Eqs. 15–17); validate these and finiteness before calling it.

The function is unit-agnostic. For solar irradiance use the same unit, for example W/m², in all three vectors. Do not silently mix power, energy, irradiance, normalized signals or time resolutions. Record the unit and forecast horizon with the result. Apply any day/night mask identically to observations and both forecasts.

## Equations and outputs

For k = 1, 2, 3, `E_k = mean(abs(y - forecast).^k)^(1/k)`. The output structure contains `E1_model` through `E3_model`, `E1_ref` through `E3_ref`, `NICE1` through `NICE3`, `NICE_Sigma`, and `n`.

`NICEk = Ek_model / Ek_persistence`; `NICE_Sigma = sum(weights .* [NICE1 NICE2 NICE3])`, with nonnegative weights summing to one for the paper's composite interpretation. Error norms retain the input unit; NICE values and the weights are dimensionless; output `n` is the count of retained samples, not the forecast horizon. Zero means perfect prediction and one means the same error as simple persistence. The paper's useful-model scale is [0,1], assuming the model is at least as good as persistence. Section 2.1 explicitly acknowledges worse models: a raw ratio above one signals this failure and must not be silently clipped. It is not a universal deployment recommendation.

## Missing data and limits

The implementation excludes a sample when any of `y`, `y_hat`, `y_ref` is NaN. It does not reject infinity, complex values or an empty retained set. Validate those cases in the caller. A zero reference error makes the ratio undefined or infinite; report that explicitly, never replace the denominator silently. The implementation does not estimate confidence intervals, prove significance, select a forecast horizon or prevent evaluation leakage. Compare models on the same retained samples.

The empirical function uses measured reference errors. `nice_theory_persistence.m` is a separate theoretical calculation with additional assumptions. See the README note about the Gaussian third absolute moment in the author manuscript: do not transfer a theoretical constant issue to the empirical ratios, and do not claim the publisher PDF has been checked if only the author version was read. Original `paper_scripts/` remain unchanged.

## Minimal example

From the repository directory, in MATLAB:

```matlab
y_full = [400 500 450 600 550];
h = 1; % one sample ahead; hourly sampling in this teaching example
y_target = y_full(1+h:end);
y_persistence = y_full(1:end-h);
m = nice_metrics(y_target, [460 480 560 580], y_persistence);
disp([m.NICE1 m.NICE2 m.NICE3 m.NICE_Sigma m.n]);
```

These are synthetic teaching irradiances in W/m², not field validation data. Expected values, rounded: 0.400, 0.365, 0.340, 0.368, 4. An agent must report whether it actually executed MATLAB; reading these values is not execution.

## Attribution

Article: C. Voyant, M. Despotovic, L. Garcia-Gutierrez, R. Amaro e Silva, P. Lauret, T. Soubdhan and N. Bailek, *NICEᵏ metrics: Unified and multidimensional framework for evaluating deterministic solar forecasting accuracy*, Sustainable Energy Technologies and Assessments 83 (2025), 104588. https://doi.org/10.1016/j.seta.2025.104588

Software concept DOI: https://doi.org/10.5281/zenodo.23562549 ; archived v1.0.0: https://doi.org/10.5281/zenodo.23562550 . Documentation added after that release is not automatically part of its archive. Cite the exact software commit or release used. License: MIT.
