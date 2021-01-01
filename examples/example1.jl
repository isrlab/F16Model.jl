using F16Model

d2r = pi/180;
npos = 0;
epos = 0;
alt = 10000;
phi = 0;
theta = 0;
psi = 0;
Vt = 300;
alp = 15*d2r;
bet = 0;
p = 0;
q = 0;
r = 0;

x0 = [npos,epos,alt,phi,theta,psi,Vt,alp,bet,p,q,r];
u0 = [9000,0,0,0,0];

# Example 1: Determine xdot
# -------------------------
# There are two ways. See below.

xdot1 = zeros(12);
F16Model.Dynamics!(xdot1,x0,u0); # This is inplace implementation -- use this with DifferentialEquations package.

xdot2 = F16Model.Dynamics(x0,u0); # Use this for linearization of dynamics, etc.

# Example 2: Least restricted trim
# --------------------------------
# Trim the aircraft for SteadyLevel at h0,V0
# xbar = trim state
# ubar = trim control
# status = status of the optimization, status = 0 means optimization found solution and (xbar, ubar) defines a valid trim  point/
# prob = data structure from IpOpt.

h0 = 10000; # ft
Vt0 = 500;   # ft/s
xbar, ubar, status, prob = F16Model.Trim(h0,Vt0); # Default is steady-level

# Example 3: Restricted trim
# ---------------------------
# Trim for steady-level and coordinated turn can be achieved by calling it with optional argumenta
# γ = flight path angle
# ψdot = turn rate
# Defaults for these are zero, as shown in example 2.
# For example, to achieve trim at 5 degrees without turning, we can call Trim(...) as
xbar, ubar, status, prob = F16Model.Trim(h0, Vt0, γ=5*pi/180, ψdot=0);

# Example 4: Most restricted trim
# -------------------------------
# Values for ϕ, ψ, θ, α, β, p, q, r are defined as pair (initial_guess, isFixed).
# The nonlinear optimization uses initial_guess to warm start the iterations. 
# If isFixed = 1, the optimization fixes the value of the state to initial_guess.
# For example, for the most restricted steady-level flight, we can call Trim(...) as follows. Note α and θ are left as free, since they must satisfy γ = θ - α. With γ = 0, we get θ = α.
xbar, ubar, status, prob = F16Model.Trim(h0, Vt0, γ=0, ψdot=0, ϕ=(0,1), ψ=(0,1), β=(0,1), p=(0,1), q=(0,1), r=(0,1)); 

# Example 5: Linearization about a given trim (xbar,ubar)
# -------------------------------------------------------
A, B = F16Model.Linearize(xbar,ubar);


