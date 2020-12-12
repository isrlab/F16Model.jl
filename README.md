# F16Model
Nonlinear model of F16 flight dynamics.

See src/testF16Model.jl for usage. 
The dynamics and linearization is functional and can be used for control law design and simulations. 

## Example

``` julia
# Add package using github url: https://github.com/isrlab/F16Model
# e.g. Pkg.add("https://github.com/isrlab/F16Model);

using F16Model

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
u0 = [9000,0,0,0,0];

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

