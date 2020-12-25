# Nonlinear Model of F16 Aircraft with NASA _table
using ForwardDiff, Ipopt
include("AeroTableInterpolation.jl")

# -------- Atmosphere Model ---------
function atmos(alt::Real,vt::Real)::Tuple{Real,Real,Real}
    rho0 = 2.377e-3;
    tfac =1 - .703e-5*(alt);

    temp = 519.0*tfac;
    if (alt >= 35000.0)
       temp=390;
    end

    rho=rho0*tfac^4.14;
    mach = (vt)/sqrt(1.4*1716.3*temp);
    qbar = .5*rho*vt^2;
    ps   = 1715.0*rho*temp;

    if ps == 0
      ps = 1715;
    end

  return (mach, qbar, ps);
end

function StateAndControlBounds()::Matrix
  d2r = pi/180;  
  bounds = [5000       40000;  # h ft
            -90*d2r    90*d2r;     # phi rad
            -90*d2r    90*d2r;     # theta rad
            -90*d2r    90*d2r;     # psi rad
             300       900;    # Vt ft/s
             -20*d2r   45*d2r;     # Alpha rad
             -30*d2r   30*d2r;     # Beta rad
            -Inf       Inf;    # p rad/s
            -Inf       Inf;    # q rad/s
            -Inf       Inf;    # r rad/s
            1000       19000;  # Thrust lbs
             -25       25;     # dele degrees
             -21.5     21.5;   # dail degrees
             -30       30;     # drud degrees
               0       25;     # dlef degrees
            ];

  return bounds;
end            


# -------- Nonlinear Dynamic of F16 model ---------
function Dynamics!(xdot::Vector, x::Vector, u::Vector)
    # Define Constants for the Aircraft Model
    g = 32.17;          # gravity, ft/s^2 
    m = 636.94;         # mass, slugs 
    B = 30.0;           # span, ft 
    S = 300.0;          # planform area, ft^2 
    cbar = 11.32;       # mean _ chord, ft 
    xcgr = 0.35;        # reference center of gravity as a fraction of cbar 
    xcg  = 0.30;        # center of gravity as a fraction of cbar. 
    Heng = 0.0;         # turbine momentum along roll axis. 

    # NasaData -- translated via eq. 2.4-6 on pg 80 of Stevens and Lewis
    Jy  = 55814.0;       # slug-ft^2  
    Jxz = 982.0;         # slug-ft^2      
    Jz  = 63100.0;       # slug-ft^2 
    Jx  = 9496.0;        # slug-ft^2 
    r2d = 180.0/pi;

    # Collect the states and control variables
    npos  = x[1];
    epos  = x[2];
    alt   = x[3];   # altitude 
    phi   = x[4];   # Euler angles are in rad. 
    theta = x[5];
    psi   = x[6];
  
    vvt   = x[7];     # total velocity 
    alpha = x[8]*r2d; # angle of attack in degrees 
    beta  = x[9]*r2d; # sideslip angle in degrees 
    P     = x[10];    # Roll Rate --- rolling  moment is Lbar 
    Q     = x[11];    # Pitch Rate--- pitching moment is M 
    R     = x[12];    # Yaw Rate  --- yawing   moment is N 
  
    sa    = sin(x[8]); # sin(alpha) 
    ca    = cos(x[8]); # cos(alpha) 
    sb    = sin(x[9]); # sin(beta)  
    cb    = cos(x[9]); # cos(beta)  
    tb    = tan(x[9]); # tan(beta)  
  
    st    = sin(theta);
    ct    = cos(theta);
    tt    = tan(theta);
    sphi  = sin(phi);
    cphi  = cos(phi);
    spsi  = sin(psi);
    cpsi  = cos(psi);

    # Control variables
    T     = u[1];   # Thrust 
    el    = u[2];   # Elevator setting in degrees 
    ail   = u[3];   # Ailerons mex setting in degrees 
    rud   = u[4];   # Rudder setting in degrees 
    lef   = u[5];   # Leading edge flap setting in degrees

    # Normalization and scaling
    if vvt <= 0.1
      vt = 0.1;
    else
      vt = vvt;
    end

    dail  = ail/21.5; # Aileron angle normalized against max deflection
    drud  = rud/30.0;  # Rudder normalized against max angle
    dlef  = (1 - lef/25.0);  # Leading edge flap normalized against max angle 

    # Atmospheric effects
    mach, qbar, ps = atmos(alt,vt);

    # Navigation equations
    U = vt*ca*cb;  # directional velocities.
    V = vt*sb;
    W = vt*sa*cb;

    nposd = U*(ct*cpsi) + V*(sphi*cpsi*st - cphi*spsi) + W*(cphi*st*cpsi + sphi*spsi); # npos dot
    eposd = U*(ct*spsi) + V*(sphi*spsi*st + cphi*cpsi) + W*(cphi*st*spsi - sphi*cpsi); # epos dot
    altd = U*st - V*(sphi*ct) - W*(cphi*ct); # alt dot

    # Kinematic equations
    phid = P + tt*(Q*sphi + R*cphi); # phi dot
    thd = Q*cphi - R*sphi; # theta dot
    psid = (Q*sphi + R*cphi)/ct; # psi dot

    # Aerodynamic forces and moments
    Cx = _Cx(alpha,beta,el); 
    Cy = _Cy(alpha,beta);  
    Cz = _Cz(alpha,beta,el);  
    Cm = _Cm(alpha,beta,el);   
    Cn = _Cn(alpha,beta,el);
    Cl = _Cl(alpha,beta,el);

    delta_Cx_lef, delta_Cy_lef, delta_Cz_lef,delta_Cm_lef, delta_Cn_lef, delta_Cl_lef  = Delta_lef(alpha,beta); 

    # Aerodynamic Damping
    Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp  = _Damping(alpha);
    delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef, delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef, delta_Cnp_lef  = _Damping_lef(alpha);

    # Rudder Influence
    delta_Cy_r30, delta_Cn_r30, delta_Cl_r30 = _RudderInfluence(alpha,beta);

    # Aileron Influence
    delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef, delta_Cl_a20, delta_Cl_a20_lef = _AileronInfluence(alpha,beta);

    # Other Coefficients
    delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el = _OtherCoefficients(alpha,el);
    delta_Cm_ds = 0; # Ignore deep-stall effects.

    # ========= Dynamics -- Euler's first and second law ==============

    # Compute total Cx 
    dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);
    Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

    # Compute total Cz
    dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);
    Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

    # Compute total Cm
    dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);    
    Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;
    
    # Compute total Cy
    dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;
    dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);
    dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);
    Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;
    
    # Compute total Cn
    dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;
    dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);
    dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);
    Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;
    
    # Compute total Cl    
    dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;
    dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);
    dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);
    Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;
    
    # Compute Udot,Vdot, Wdot,(as on NASA report p36)
    Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;
    Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;
    Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;
    
    # vt_dot equation 
    vd = (U*Udot + V*Vdot + W*Wdot)/vt;
    
    # alpha_dot equation
    αd = (U*Wdot - W*Udot)/(U*U + W*W);
    
    # beta_dot equation
    βd = (Vdot*vt - V*vd)/(vt*vt*cb);
    
    # Equations for Pdot, Qdot, and Rdot 
    L_tot = Cl_tot*qbar*S*B;       # Get moments from coefficients
    M_tot = Cm_tot*qbar*S*cbar;
    N_tot = Cn_tot*qbar*S*B;
    denom = Jx*Jz - Jxz*Jxz;
    
    # Pdot
    pd =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;
    
    # Qdot
    qd = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;
    
    # Rdot
    rd = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

    xdot[1] = nposd;
    xdot[2] = eposd;
    xdot[3] = altd;
    xdot[4] = phid;
    xdot[5] = thd;
    xdot[6] = psid;
    xdot[7] = vd;
    xdot[8] = αd;
    xdot[9] = βd;
    xdot[10] = pd;
    xdot[11] = qd;
    xdot[12] = rd;
end

# Redefine the function for gradient calculation.
function Dynamics(x::Vector, u::Vector)::Vector
  # Define Constants for the Aircraft Model
  g = 32.17;          # gravity, ft/s^2 
  m = 636.94;         # mass, slugs 
  B = 30.0;           # span, ft 
  S = 300.0;          # planform area, ft^2 
  cbar = 11.32;       # mean _ chord, ft 
  xcgr = 0.35;        # reference center of gravity as a fraction of cbar 
  xcg  = 0.30;        # center of gravity as a fraction of cbar. 
  Heng = 0.0;         # turbine momentum along roll axis. 

  # NasaData -- translated via eq. 2.4-6 on pg 80 of Stevens and Lewis
  Jy  = 55814.0;       # slug-ft^2  
  Jxz = 982.0;         # slug-ft^2      
  Jz  = 63100.0;       # slug-ft^2 
  Jx  = 9496.0;        # slug-ft^2 
  r2d = 180.0/pi;

  # Collect the states and control variables
  npos  = x[1];
  epos  = x[2];
  alt   = x[3];    # altitude 
  phi   = x[4];   # Euler angles are in rad. 
  theta = x[5];
  psi   = x[6];

  vvt   = x[7];     # total velocity 
  alpha = x[8]*r2d; # angle of attack in degrees for _ tables
  beta  = x[9]*r2d; # sideslip angle in degrees for _ tables
  P     = x[10];    # Roll Rate --- rolling  moment is Lbar 
  Q     = x[11];    # Pitch Rate--- pitching moment is M 
  R     = x[12];    # Yaw Rate  --- yawing   moment is N 

  sa    = sin(x[8]); # sin(alpha) 
  ca    = cos(x[8]); # cos(alpha) 
  sb    = sin(x[9]); # sin(beta)  
  cb    = cos(x[9]); # cos(beta)  
  tb    = tan(x[9]); # tan(beta)  

  st    = sin(theta);
  ct    = cos(theta);
  tt    = tan(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);
  spsi  = sin(psi);
  cpsi  = cos(psi);

  # Control variables
  T     = u[1];   # Thrust 
  el    = u[2];   # Elevator setting in degrees 
  ail   = u[3];   # Ailerons mex setting in degrees 
  rud   = u[4];   # Rudder setting in degrees 
  lef   = u[5];   # Leading edge flap setting in degrees

  # Normalization and scaling
  if vvt <= 0.1
    vt = 0.1;
  else
    vt = vvt;
  end

  dail  = ail/21.5; # Aileron angle normalized against max deflection
  drud  = rud/30.0;  # Rudder normalized against max angle
  dlef  = (1 - lef/25.0);  # Leading edge flap normalized against max angle 

  # Atmospheric effects
  mach, qbar, ps = atmos(alt,vt);

  # Navigation equations
  U = vt*ca*cb;  # directional velocities.
  V = vt*sb;
  W = vt*sa*cb;

  nposd = U*(ct*cpsi) + V*(sphi*cpsi*st - cphi*spsi) + W*(cphi*st*cpsi + sphi*spsi); # npos dot
  eposd = U*(ct*spsi) + V*(sphi*spsi*st + cphi*cpsi) + W*(cphi*st*spsi - sphi*cpsi); # epos dot
  altd = U*st - V*(sphi*ct) - W*(cphi*ct); # alt dot

  # Kinematic equations
  phid = P + tt*(Q*sphi + R*cphi); # phi dot
  thd = Q*cphi - R*sphi; # theta dot
  psid = (Q*sphi + R*cphi)/ct; # psi dot

  # Aerodynamic forces and moments
  Cx = _Cx(alpha,beta,el); 
  Cy = _Cy(alpha,beta);  
  Cz = _Cz(alpha,beta,el);  
  Cm = _Cm(alpha,beta,el);   
  Cn = _Cn(alpha,beta,el);
  Cl = _Cl(alpha,beta,el);

  delta_Cx_lef, delta_Cy_lef, delta_Cz_lef,delta_Cm_lef, delta_Cn_lef, delta_Cl_lef  = Delta_lef(alpha,beta); 

  # Aerodynamic Damping
  Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp  = _Damping(alpha);
  delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef, delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef, delta_Cnp_lef  = _Damping_lef(alpha);

  # Rudder Influence
  delta_Cy_r30, delta_Cn_r30, delta_Cl_r30 = _RudderInfluence(alpha,beta);

  # Aileron Influence
  delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef, delta_Cl_a20, delta_Cl_a20_lef = _AileronInfluence(alpha,beta);

  # Other Coefficients
  delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el = _OtherCoefficients(alpha,el);
  delta_Cm_ds = 0; # Ignore deep-stall effects.

  # ========= Dynamics -- Euler's first and second law ==============

  # Compute total Cx 
  dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);
  Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

  # Compute total Cz
  dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);
  Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

  # Compute total Cm
  dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);    
  Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;
  
  # Compute total Cy
  dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;
  dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);
  dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);
  Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;
  
  # Compute total Cn
  dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;
  dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);
  dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);
  Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;
  
  # Compute total Cl    
  dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;
  dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);
  dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);
  Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;
  
  # Compute Udot,Vdot, Wdot,(as on NASA report p36)
  Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;
  Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;
  Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;
  
  # vt_dot equation 
  vd = (U*Udot + V*Vdot + W*Wdot)/vt;
  
  # alpha_dot equation
  αd = (U*Wdot - W*Udot)/(U*U + W*W);
  
  # beta_dot equation
  βd = (Vdot*vt - V*vd)/(vt*vt*cb);
  
  # Equations for Pdot, Qdot, and Rdot 
  L_tot = Cl_tot*qbar*S*B;       # Get moments from coefficients
  M_tot = Cm_tot*qbar*S*cbar;
  N_tot = Cn_tot*qbar*S*B;
  denom = Jx*Jz - Jxz*Jxz;
  
  # Pdot
  pd =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;
  
  # Qdot
  qd = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;
  
  # Rdot
  rd = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

  xdot = [nposd, eposd, altd, phid, thd, psid, vd, αd, βd, pd, qd, rd];

  return xdot;
end

function Linearize(x0::Vector,u0::Vector)::Tuple{Matrix,Matrix}
  f1 = x -> Dynamics(x,u0);
  f2 = u -> Dynamics(x0,u);
  A = ForwardDiff.jacobian(f1,x0);
  B = ForwardDiff.jacobian(f2,u0);
  return (A,B)
end

#= Four standard trim routines
1. SteadyLevel -- given (h,Vt), α, θ, and all control variables are free, all xdots are zero.
2. SteadyPullUp -- given (h,Vt,θdot), ... 
3. SteadyRoll -- given (h,Vt,ϕdot), ...
4. SteadyTurning -- given (h,Vt,ψdot), ...
=#
function Trim(h,V,TrimType;phidot=0,thetadot=0,psidot=0)
  d2r = pi/180;
  h0 = h;
  Vt0 = V;
  p0,q0,r0 = 0,0,0;
  xdot_ref = zeros(10);

  if TrimType == :SteadyLevel
    phi0, theta0,psi0 = 0,0,0;
    alpha0,beta0 = 0,0;
    T0, dele0, dail0, drud0, dlef0 = 9000,0,0,0,0;
    ix = [1, 1,0,1, 1,0,1, 1,1,1]; # 1 => Trim values are fixed in the optimization, 0 => Trim values are optimization variables.
  elseif TrimType == :SteadyTurning
    error(":SteadyTurning not implemented.")
  elseif TrimType == :SteadyPullUp
    error(":SteadyPullUp not implemented.")
  elseif TrimType == :SteadyRoll
    error(":SteadyRoll not implemented.")
  else
    error("Invalid trim type specified.")
  end

  x0 = [h0, phi0, theta0, psi0, Vt0, alpha0, beta0, p0, q0, r0];
  u0 = [T0,dele0,dail0,drud0,dlef0];
  ii = findall(x->x==0, ix);

  # Extract states and control from optimization variables
  function getXU(z)
    z0 = x0
    z0[ii] .= 0;

    Mz = zeros(length(x0),length(ii));
    for (i,ind) in enumerate(ii)
        Mz[ind,i] = 1;
    end

    v = z0 + Mz*z[1:length(ii)];

    X = [0;0;v];
    U = z[(length(ii)+1):end];
    return (X,U)
  end

  # Nonlinear constraints from dynamics and trim types
  function _nonlinear_constraints(z)
    var = z[1:end-1];
    X,U = getXU(var);
    xdot = Dynamics(X, U)[3:12];
    e = xdot - xdot_ref;
    g = e'*e - z[end];
    return g
  end

# =================== IpOpt Specific Functions =====================
  function eval_f(z) # Cost function
    return z[end];
  end
  
  function eval_g(z, g) # Nonlinear inequality constraint
    # Bad: g    = zeros(2)  # Allocates new array
    # OK:  g[:] = zeros(2)  # Modifies 'in place'
    g[:] .= _nonlinear_constraints(z);
  end
  
  function eval_grad_f(z, gradf) # Gradient of cost function.
    # Bad: grad_f    = zeros(4)  # Allocates new array
    # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
    f = x-> ForwardDiff.gradient(eval_f,x);
    gradf[:] = f(z);
  end
  
  function eval_jac_g(z, mode, rows, cols, values) # Jacobian of nonlinear inequality constraint.
    nvar = length(z);
    ng = 1;

    if mode == :Structure
        index = 1;
        for i=1:ng
            for j=1:nvar
                rows[index] = i; 
                cols[index] = j;
                index += 1;
            end
        end
    else
        M = x->ForwardDiff.gradient(_nonlinear_constraints,x);
        values[:] = M(z);
    end
  end
  
  function eval_h(z, mode, rows, cols, obj_factor, lambda, values) # Hessian of cost.
    nvar = length(z);
    if mode == :Structure
      # Symmetric matrix, fill the lower left triangle only
      idx = 1
      for row = 1:nvar
        for col = 1:row
          rows[idx] = row
          cols[idx] = col
          idx += 1
        end
      end
    else
      H = obj_factor*zeros(nvar,nvar);
      f(z) = _nonlinear_constraints(z); # Only one constraint
      ∇2f = x->ForwardDiff.hessian(f,x);
      H = H + lambda[1]*∇2f(z);

      # Copy them into return variable -- only lower triangle.
      idx = 1;
      for row = 1:nvar
        for col = 1:row
          values[idx] = H[row,col];
          idx += 1;
        end
      end
    end
  end

  # Call IpOpt to solve the constraint nonlinear optmization problem
  bounds = StateAndControlBounds();
  bx = bounds[ii,:]; # State bounds
  bu = bounds[11:15,:]; # Control bounds

  n = length(ii) + 5 + 1;
  x_L = [bx[:,1];bu[:,1];0];
  x_U = [bx[:,2];bu[:,2];Inf];
  m = 1; # Only one nonlinear constraint ... true for all trim types
  g_L = [-Inf];
  g_U = [0.0];

  nele_jac  = m*n; 
  nele_hess = Integer(n*(n+1)/2);
  prob = createProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);

  prob.x = [x0[ii];u0;0]; # Initial guess.
  status = solveProblem(prob); 
  xbar,ubar = getXU(prob.x);
  return (xbar, ubar, status, prob.obj_val);
end

