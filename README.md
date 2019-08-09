# Lanczos-Filter
Low pass filtering a time-series by applying a weighted running mean over the time dimension.<br>
Details of Lanczos Filter could be found [here](https://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281979%29018%3C1016%3ALFIOAT%3E2.0.CO%3B2)

## Core Code:
**lanczosbp.py:**<br>Use two low pass lanczos filters to get 3 to 10 days bandpass 850 hPa vorticity, the variance of which could be thought of the pre-TC synoptic disturbance, i.e. TC seed index.

## Data:
EAR5 hourly 850hPa vorticity.
