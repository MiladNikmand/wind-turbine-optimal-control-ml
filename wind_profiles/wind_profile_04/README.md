# wind profile 04

Wind prediction benchmark. Generated automatically -- do not edit by hand.

## Run settings

| Parameter | Value |
|---|---|
| History window | 15.00 s |
| Forecast horizon | 1.50 s |
| Stride | 1.50 s |
| Segments | 190 |
| Wind range | 10.00 to 15.00 m/s |
| Selection criterion | RMSE |

## Model ranking

| Rank | Model | RMSE | nRMSE | MAE | MAPE % | sMAPE % | R2 | DM p | Time (s) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `rbf` | 0.0511 | 0.0102 | 0.0451 | 0.37 | 0.37 | 0.7807 | 1.0000 | 0.01 |
| 2 | `rbf_arima` | 0.0511 | 0.0102 | 0.0451 | 0.37 | 0.37 | 0.7807 | NaN | 0.08 |
| 3 | `arima` | 0.1479 | 0.0296 | 0.1114 | 0.90 | 0.90 | -2.6181 | 0.0000 | 0.29 |

## Per-model forecasts

### RBF (rank 1, RMSE 0.0511)

![rbf](prediction/web/figures/rbf_full.png)

### RBF_ARIMA (rank 2, RMSE 0.0511)

![rbf_arima](prediction/web/figures/rbf_arima_full.png)

### ARIMA (rank 3, RMSE 0.1479)

![arima](prediction/web/figures/arima_full.png)

## Artifacts

- `report/prediction_suite_report.txt` -- full written report
- `report/prediction_suite_summary.json` -- machine-readable summary
- `report/*.csv` -- metrics as CSV (GitHub renders these as tables)
- `summary_table.tex` -- LaTeX table for the thesis
- `prediction/figures/`, `prediction/animation/` -- full-resolution assets
- `prediction/web/` -- downscaled assets used by report.html
- `prediction_suite_results.mat` -- raw numbers behind all of the above
