# rsv_model_fitting
Code for fitting SEIRS transmission models for respiratory syncytial virus (RSV) to incidence data.

## Authors
Code authors: Alexandra B Hogan, Roslyn I Hickson, Geoffry N Mercer, Hannah C Moore, Kathryn Glass (2016).

## About the model
The model here is an ordinary differential equation (ODE) transmission model for RSV transmission for two age classes. This model is published in:  

Hogan AB, Glass K, Moore HC, Anderssen R. Exploring the dynamics of respiratory syncytial virus (RSV) transmission in children. *Theor. Popul. Biol.* 2016;110:78–85. Available from: http://dx.doi.org/10.1016/j.tpb.2016.04.003.

The parameters for ageing, birth rates, latency and infectious period duration, are all RSV-specific and specified at two different sections in the code. The timestep is currently in weeks but this can be changed.

## About the fitting  
The fitting is implemented using the Nelder Mead algorithm, via the fminsearch function in MATLAB. You will need to start with reasonable parameter guesses for the fitting to work. There is an option in the code to either run the model with the parameter guesses, or run the fitting routine. The fitted parameters are the transmission coefficients for the two age groups; the amplitude of seasonal forcing; the phase shift (used to obtain horizontal alignment); and the duration of immunity.

## Synthetic data
The fitting routine is demonstrated using synthetic data. Note that this 'data' was generated from the model itself, with a negative binomial distribution used to add in noise, and is therefore only for illustrative purposes.
