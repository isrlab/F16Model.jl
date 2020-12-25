# F16Model

This is a Julia package for a nonlinear model of the F16 aircraft. The aerodynamics included in this model come from the NASA Technical Report 1538, *Simulator Study of Stall/Post-Stall Characteristics of a Fighter Airplane with Relaxed Longitudinal Static Stability*, by Nguyen, Ogburn, Gilbert, Kibler, Brown, and Deal, Dec 1979. The flight dynamics are based on *Aircraft Control and Simulations*, by Brian Stevens and Frank Lewis, Wiley Inter-Science, New York, 1992. 

This Julia package is inspired by the [MATLAB/Simulink package](https://dept.aem.umn.edu/~balas/darpa_sec/SEC.Software.html), and currently has the following features:

1. Nonlinear dynamics for simulations.
2. Linearization about a trim point (x0,u0).

More features will be added as we continue to develop this package.

## Model Details 
Detailed information is available [here](https://dept.aem.umn.edu/~balas/darpa_sec/software/F16Manual.pdf).

The 12 states of the system are as follows:

1. N: North position in ft
1. E: East position in ft
1. h: Altitude in ft, min: 5000 ft, max: 40000 ft
1. phi: Roll angle in rad
1. theta: Pitch angle in rad
1. psi: Yaw angle in rad
1. Vt: Magnitude of total velocity in ft/s, min: 300 ft/s, max: 900 ft/s
1. alpha: Angle of attack in rad, min: -20 deg, max: 45 deg
1. beta: Side slip angle in rad, min: -30 deg, max: 30 deg
1. p: Roll rate in rad/s
1. q: Pitch rate in rad/s
1. r: Yaw rate in rad/s

The 5 control variables are:

1. T: Thrust in lbs, min: 1000, max: 19000 
1. dele: Elevator angle in deg, min:-25, max: 25
1. dail: Aileron angle in deg, min:-21.5, max: 21.5
1. drud: Rudder angle in deg, min: -30, max: 30
1. dlef: Leading edge flap in deg, min: 0, max: 25

Actuator models are defined as:

1. T: max |rate|: 10,000 lbs/s 
1. dele: max |rate|: 60 deg/s
1. dail: max |rate|: 80 deg/s
1. drud: max |rate|: 120 deg/s
1. dlef: max |rate|: 25 deg/s

The nonlinear model of the aircraft does not include actuator dynamics. 
The actuator dynamics need to be modeled as LTI systems and added to the system.
For example, for dele the low pass filter 1/(s/60+1) would model the actuator dynamics.

## Installation
Add package using GitHub url as shown below.

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
alt = 10000; # should be in between 5000 ft and 100000 ft
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
F16Model.Dynamics!(xdot1,x0,u0); # Does not implement actuator dynamics.

# Evaluate xdot -- returns vector. Use this for linearization of dynamics, etc.
xdot2 = F16Model.Dynamics(x0,u0); #  Does not implement actuator dynamics.

# Linearize about some trim point (x0,u0)
A, B = F16Model.Linearize(x0,u0);
```

## Trim Functions

The aircraft model can be trimmed at the following configurations:

1. **Steady-level Flight:** phi = 0, phidot = 0, thetadot = 0, psidot = 0
2. **Steady Turning Flight:** phidot = 0, thetadot = 0, psidot = given
3. **Steady Pull Up:** phi = 0, phidot = 0, psidot = 0, thetadot = given
4. **Steady Roll:** thetadot = 0, psidot = 0, phidot = given

For all flight conditions: pdot, qdot, udot, Vtdot, alphadot, and betadot are all zero. Derivatives of states N and E are ignored in the trim calculations.

The trim conditions are determined by solving a constrained nonlinear optimization problem.
min (xdot-xdot_ref)*(xdot-xdot_ref) subject to state and control constraints. 

The optimization returns the tuple (xbar, ubar, status, objVal). 
If status == 0, the nonlinear optimization was succesful.
If objVal is small then (xbar,ubar) are valid trim state and control values.

For example:

```julia
# Trim the aircrat

h0 = 10000; # Trim at this altitude
Vt0 = 500;  # Trim at this velocity

# Trim for steady level flight ay height h0 and velocity Vt0
# ----------------------------------------------------------
xbar, ubar, status, objVal = F16Model.Trim(h0,Vt0,:SteadyLevel); 