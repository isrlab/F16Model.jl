# F16Model

This is a Julia package for a nonlinear model of the F16 aircraft. The aerodynamics included in this model come from the NASA Technical Report 1538, *Simulator Study of Stall/Post-Stall Characteristics of a Fighter Airplane with Relaxed Longitudinal Static Stability*, by Nguyen, Ogburn, Gilbert, Kibler, Brown, and Deal, Dec 1979. The flight dynamics are based on *Aircraft Control and Simulations*, by Brian Stevens and Frank Lewis, Wiley Inter-Science, New York, 1992. This Julia package aims to replicate the MATLAB/Simulink package developed by Prof. Gary J. Balas with his students from the Aerospace Engineering & Mechanics department of the University of Minnesota and his collaborators. We also hope to add more features to the Julia package.

Currently, the package supports the following features:

1. Nonlinear dynamics for simulations.
2. Lineariation about a trim point (x0,u0).

More features will be added as we continue to develop this package.

## Installation
Add package using github url: https://github.com/isrlab/F16Model, e.g.

``` julia
Pkg.add("https://github.com/isrlab/F16Model");
```

## Example

``` julia
using F16Model

# Define state vector
# -------------------
d2r = pi/180;
npos = 0;
epos = 0;
alt = 10000;
phi = 0;
theta = 0;
psi = 0;
Vt = 300;
alp = 0;
bet = 0;
p = 0;
q = 0;
r = 0;

x0 = [npos,epos,alt,phi,theta,psi,Vt,alp,bet,p,q,r];

# Define control vector
# ---------------------
T = 9000; # Thrust lbs
dele = 0; # deg elevator angle
dail = 0; # deg aileron angle
drud = 0; # deg rudder angle
dlef = 0; # deg leading edge flap angle
u0 = [T,dele,dail,drud,dlef];

# Evaluate xdot -- inplace implementation -- use this with DifferentialEquations package.
xdot1 = zeros(12);
F16Model.Dynamics!(xdot1,x0,u0);

# Evaluate xdot -- returns vector
xdot2 = F16Model.Dynamics(x0,u0); # Use this for linearization of dynamics, etc.

# Linerize about some trim point (x0,u0)
A, B = F16Model.Linearize(x0,u0);
```

## To do:
1. Trim functions
2. Aero-table plotting

