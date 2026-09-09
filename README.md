# wind-turbine-optimal-control-ml

Master's thesis project on wind turbine MPPT control. Wind speed is first predicted using machine learning, then a proposed DDP-HJB solver computes the optimal control offline over that predicted horizon — so the optimal solution is ready and applied in real time, enabling near-optimal performance without heavy onboard computation.

**[View the full results report →](https://miladnikmand.github.io/wind-turbine-optimal-control-ml/report.html)**

---

## How it works

The pipeline runs in two stages over the same segmentation geometry.

**1. Segmented wind prediction.** A wind profile is sliced into overlapping windows: `history_sec` of training data, then a `predict_sec` forecast horizon, advancing by `stride_sec`. Ten forecasting models each predict every segment. They are scored per segment on RMSE, MAE, MAPE, sMAPE, nRMSE and R², ranked, and compared against the best model with a Diebold-Mariano test. A per-segment winner is selected and the winners are stitched into one continuous predicted series.

**2. DDP-HJB optimal control.** The controller reads those cached forecasts and solves the optimal control problem over each predicted segment, carrying state across segment boundaries. The turbine is a four-state single-mass model, `x = [Wr; lambda; Cp; Tg]`, driven by a single scalar control input acting on the generator-torque dynamics. Because the solve happens over a forecast rather than measured wind, the control is ready before the wind arrives.

The forecast source is selectable, which is what makes the cost of prediction error measurable:

| Mode | Meaning |
|---|---|
| `truth` | Perfect foresight — the upper bound |
| `bestof` | The per-segment winning model |
| `model` | One named predictor, e.g. `rbf` |

The gap between `truth` and the others is the price of forecast error.

Both stages report a real-time budget: mean solve time per segment against the forecast horizon, so you can see whether the method fits inside its own deadline.

---

## Requirements

**MATLAB R2020a or later.** The floor is `exportgraphics`, used unguarded in `+report`. No other function in the project requires anything newer.

Eight of the ten predictors run on **base MATLAB with no toolboxes**. The rest degrade gracefully — `predictors.available()` checks what your install can run and `run_prediction` drops the unavailable ones rather than failing partway through a batch.

| Toolbox | Needed by |
|---|---|
| Statistics and Machine Learning | `rbf`, `rbf_arima`, `svr`, `gp` |
| Econometrics | `arima`, `rbf_arima` |
| Deep Learning **(R2023b or earlier)** | `bilstm` |

`bilstm` calls `trainNetwork`, which MathWorks removed in R2024a. On a newer release it is skipped automatically; use `gru`, `lstm` or `tcn` instead.

---

## Quick start

Everything is driven from the project root and there are no interactive prompts.

```matlab
>> cd path/to/wind-turbine-optimal-control-ml
>> check_packages        % is every file where MATLAB expects it?
>> run_all               % prediction + control, one profile, smoke-run sized
```

`run_all` ships configured for a short run: one profile, three predictors, five control segments. Widen it by editing the CONFIG block once you've confirmed the report looks right.

To run the stages separately:

```matlab
>> profile_id = 'wind03';
>> run_prediction        % writes prediction_suite_results.mat and all figures
>> run_control           % solves DDP-HJB over the cached forecasts
```

Then open `report.html`.

**Run from the project root.** MATLAB resolves `+package` folders and the `wind_profiles/` output path relative to the current folder.

**`clear` first if running standalone.** `run_prediction` and `run_control` guard every setting with `if ~exist(...,'var')` so `run_all` can drive them. That also means stale variables from an earlier session silently override the defaults.

---

## Layout

```
+predictors/     ten forecasting models + registry, availability gate,
                 segment geometry, stitching, and the forecast cache
                 that feeds the control stage
+ddp/            DDP-HJB solver: forward pass, backward pass, RK4
                 integrator, second-order tensors, reference
                 trajectories, and all tuning constants in params.m
+windlib/        wind profile registry, loader, resampler, feasibility
                 checker, and the source CSVs in +windlib/data/
+report/         figure export, GIF animation, statistics plots, and the
                 JSON/JS bundle report.html reads

run_all.m        drives the whole pipeline unattended, skipping
                 completed work
run_prediction.m prediction benchmark for one profile
run_control.m    DDP-HJB control for one profile and forecast mode
report.html      the results viewer

wind_profiles/   generated outputs, one folder per profile
```

Supporting scripts: `check_packages`, `check_report` and `diagnose_wind` verify an install; `test_wind`, `test_predictors`, `test_ddp` and `test_report` are per-package smoke tests; `Bench_rbf` profiles where RBF spends its time.

### Adding a predictor

Drop a `.m` file in `+predictors/` matching the shared contract, then add one line to the table in `registry.m`. Nothing else changes.

```matlab
[future_pred, metrics, segment_time_sec] = fn( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
```

---

## Wind profiles

Seven profiles, referred to by id or by number.

| id | Description | Duration |
|---|---|---|
| `wind01` | Original, 12–24 m/s | 300 s |
| `wind02` | Small wind data | 60 s |
| `wind03` | Third small wind data | 60 s |
| `wind04` | `wind01` rescaled to 10–15 m/s | 300 s |
| `wind05` | Synthetic, `12 + sin(0.7t)`, smoothed | 100 s |
| `wind06` | IEEE paper wind | 10 s |
| `wind07` | Renewable-Energy 2016 | 100 s |

`wind06` is only 10 s long and yields **no segments** at the default `history_sec = 15`. It needs `history_sec <= 8.5` at a 1.5 s horizon. Run `windlib.feasibility(step_size, history_sec, predict_sec, stride_sec)` to see the segment count for every profile before committing to a long sweep.

---

## Results report

`report.html` is a single self-contained viewer that reads `wind_profiles/report_data.js`, regenerated by `report.build_data` after every stage. It has to sit next to the `wind_profiles/` folder.

GitHub renders HTML as source, so use the Pages link at the top of this file rather than opening the file in the repository browser. Opening `report.html` locally also works, though the JSON fallback path is blocked by browser CORS rules over `file://` — the `.js` bundle is there for exactly that reason. MathJax is loaded from a CDN, so equations need an internet connection either way.

---

## Configuration notes

**The two entry points default to different `step_size`.** `run_prediction` standalone uses `0.01`; `run_all` sets `0.1`. Under `run_all` the controller step `dt_controller = 0.01` no longer matches, so the segment resample becomes a 10× upsample rather than an identity. Set both deliberately.

**`skip_done` skips on file existence, not on settings match.** If you change `history_sec`, `predict_sec` or `stride_sec` and re-run with `skip_done = true`, prediction is skipped because the `.mat` exists, and the control stage then refuses to start because the segment indices no longer line up. Delete that profile's `prediction_suite_results.mat` after changing any window setting.

**Segment figures are off by default.** `export_segment_figures` writes roughly four PNGs per model per segment — about 7,600 files for `wind01`. The `.gitignore` excludes that tier.

---

## Known issues

This project is a faithful refactor of the original thesis code. Several defects in the source were reproduced deliberately rather than corrected, so that the numbers remain reproducible. They are documented here and in `+ddp/params.m` rather than silently fixed.

**Control cost uses the rotor radius.** The backward pass computes `l_u = R * U` and `l_uu = R`, where `R` is the rotor radius (4.4), not the control weight `Rmat`. `Rmat` is used only in the reported `true_cost`, so the cost the gains are derived from disagrees with the cost that gets logged.

**Dead cost-weight assignment.** `Q` and `Qf` were assigned twice in the original; the first pair was immediately overwritten. The weights actually in use are the `Q1..Q4` set. The dead assignment is preserved as a comment.

**Row-vector `mrdivide` in the A matrix.** `A32` divides two row vectors, collapsing to a scalar. The intended form is noted inline in `+ddp/backward.m`.

**Two different Cp models.** The RK4 stages evaluate `Cp = 0.22*(116/lambda - 5)*exp(-12.5/lambda)`, while the derivative terms use `exp(0.4375/(beta^3+1) - 12.5/D)`. These disagree.

**Second-order tensors are inconsistent with the state ordering.** `ddp.second_order` reads `x(3)` as `Tg` and `x(4)` as `Cp`, the reverse of the ordering used everywhere else, and writes its tensor rows to match. Since `Tg_dot` is affine, its second derivatives must be zero, yet the shipped tensor places Cp curvature there; finite differencing confirms the pattern from the opposite direction. The tensor also contains entries that are first derivatives rather than second. Measured effect on the solve: tracking RMSE and cost each shift by **0.004%**, with `max |dWr|` of 0.0029 rad/s against an 18.4 rad/s operating point. The second-order contribution is heavily damped by `alpha`, `alpha_filter` and the `P`/`s` caps before it reaches the gains, so results are unaffected. `test_second_order.m` reproduces the diagnosis.

**Aerodynamic torque differs between the integrator and the Jacobian.** `A12` and `A13` in `+ddp/backward.m` imply `Ta ∝ Cp/lambda³`, while `+ddp/rk4.m` computes `Ta = 0.5·ρ·π·R³·v²·Cp/lambda`. This one is unresolved and sits in the first-order gains.

**Turbine constants are duplicated.** `+ddp/rk4.m` and `+ddp/second_order.m` hardcode their own copies of `R`, `Kt`, `Jt`, `Rs`, `Rho` and `Ls` instead of taking them from `ddp.params`. They agree today; editing `params.m` alone will silently desync them.

---

## License

Released under the MIT License — see [LICENSE](LICENSE). You are free to use, modify and redistribute this work, including commercially, provided the copyright notice and license text are retained.

## Citation

If you use this work, please cite it. See [CITATION.cff](CITATION.cff), or use GitHub's "Cite this repository" button.

> [Author]. *[Thesis title]*. Master's thesis, [Institution], [Year].
