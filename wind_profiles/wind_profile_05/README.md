# wind profile 05

Wind prediction benchmark. Generated automatically -- do not edit by hand.

## Run settings

| Parameter | Value |
|---|---|
| History window | 15.00 s |
| Forecast horizon | 1.50 s |
| Stride | 1.50 s |
| Segments | 56 |
| Wind range | 11.58 to 12.72 m/s |
| Selection criterion | RMSE |

## Model ranking

| Rank | Model | RMSE | nRMSE | MAE | MAPE % | sMAPE % | R2 | DM p | Time (s) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `arima` | 0.0428 | 0.0375 | 0.0322 | 0.27 | 0.27 | -9.5507 | 1.0000 | 0.23 |
| 2 | `rbf` | 0.5568 | 0.4880 | 0.3354 | 2.78 | 3.10 | -465.7415 | 0.0000 | 0.01 |
| 3 | `rbf_arima` | 0.5568 | 0.4880 | 0.3354 | 2.78 | 3.10 | -465.7415 | 0.0000 | 0.08 |

## Per-model forecasts

### ARIMA (rank 1, RMSE 0.0428)

![arima](prediction/web/figures/arima_full.png)

### RBF (rank 2, RMSE 0.5568)

![rbf](prediction/web/figures/rbf_full.png)

### RBF_ARIMA (rank 3, RMSE 0.5568)

![rbf_arima](prediction/web/figures/rbf_arima_full.png)

## Artifacts

- `report/prediction_suite_report.txt` -- full written report
- `report/prediction_suite_summary.json` -- machine-readable summary
- `report/*.csv` -- metrics as CSV (GitHub renders these as tables)
- `summary_table.tex` -- LaTeX table for the thesis
- `prediction/figures/`, `prediction/animation/` -- full-resolution assets
- `prediction/web/` -- downscaled assets used by report.html
- `prediction_suite_results.mat` -- raw numbers behind all of the above
