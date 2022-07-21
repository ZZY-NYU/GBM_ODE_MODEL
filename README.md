Code for GBM simulation
========

Files
--------
* `Fig5H+Fig5S3_corrected.xlsx` contains the data of apoptosis, CD154 and CD163 from Xin Cui's experiment. Some of the data were obviously affected by experiment inaccuracy, and was masked out. The masked data were labeled in `Fig5H+Fig5S3_digest.xlsx`

* `Fig5F+Fig5S2+Fig3G.xlsx` contains cytokine data, but was not used for manuscript. It is placed here only for code dependency reason.

* `code/PD_PDL1_Integral_Inputs_ATFL_influx_slope_us_problem.py` is the equations that determine the mathematical structure of the system.

* `code/new_fit33.py` is the script to fit the parameters.

* `code/fit_check_new_band.py` is the script to check the validation patient group with the fitted parameters and generate the figures used in the manuscript.

* `code/log.d/` is to contain log files which take a snap of the scripts every time they are run.

* `code/fig_exp.d/new_fit/` contains the fitted parameters. (The parameters file used for manuscript is included inside.)

* `code/fig_exp.d/fit_check/` is to contain the generated figures.


How to run the code
--------

### Dependencies
Please make sure you have Python 3.x installed with the listed modules available
* numpy
* matplotlib
* scipy
* openpyxl

### To check the current parameters
The fitted parameters used for manuscript is `code/fig_exp.d/new_fit/exp_20210202_210953_res.bin`. To seee the parameters, run in Python

```python
import pickle
with open('code/fig_exp.d/new_fit/exp_20210202_210953_res.bin','rb') as f:
    params = pickle.load(f)[1]
print(params)
```

### To run the simulation
To run the simualtion with current parameters, just try repeating the current `fit_check_new_band.py`. The way to use the model is mostly self-explanatory if you read the code.

### To fit the parameters yourself
If you have your own data, or you want to play with the fitting, run `new_fit33.py`. It will take a while, and if the solver converges, the result will be saved in `code/fig_exp.d/new_fit/` with current datetime in the filename.


