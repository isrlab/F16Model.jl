include("NonlinearF16Model.jl");

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

# Generate random state vector within constraints.
# xdot = zeros(12);
# alt = rand(7000:10000);
# phi = rand(-pi/3:pi/3);
# theta = rand(-pi/3:pi/3);
# psi = rand(-pi/3:pi/3);
# Vt = rand(1000.0:9000.0);
# alp = rand(-20:45)*d2r; # Rad
# bet = rand(-30:30)*d2r; # Rad
# p = rand(-20.0:20.0);
# q = rand(-20.0:20.0);
# r = rand(-20.0:20.0);

# T = rand(1000:15000.0);
# el = rand(-25.0:25.0);
# ail = rand(-20.0:20.0);
# rud = rand(-30.0:30.0);
# lef = rand(-25:25.0);

x0 = [npos,epos,alt,phi,theta,psi,Vt,alp,bet,p,q,r];
u0 = [9000,0,0,0,0];

xdot = zeros(12);

# Check speed of implementation 1
@time for i=1:10000
    Dynamics!(xdot,x0,u0);
end

# Check speed of implementation 2
@time for i=1:10000
    xdot = Dynamics(x0,u0);
end

A, B = Linearize(x0,u0);

# Todo
# ----
# 1. Write trim routine and solve it using nonlinear optimization.
# 2. Generate sensor models -- output acceleration, dynamic pressure, etc.
# 3. Disturbance models.
# 4. Plot aero tables.
