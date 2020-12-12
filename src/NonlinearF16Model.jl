# Nonlinear Model of F16 Aircraft with NASA Aerotable
using ForwardDiff
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

# -------- Nonlinear Dynamic of F16 model ---------
function Dynamics!(xdot::Vector, x::Vector, u::Vector)
    # Define Constants for the Aircraft Model
    g = 32.17;          # gravity, ft/s^2 
    m = 636.94;         # mass, slugs 
    B = 30.0;           # span, ft 
    S = 300.0;          # planform area, ft^2 
    cbar = 11.32;       # mean aero chord, ft 
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
    Cx = AeroCx(alpha,beta,el); 
    Cy = AeroCy(alpha,beta);  
    Cz = AeroCz(alpha,beta,el);  
    Cm = AeroCm(alpha,beta,el);   
    Cn = AeroCn(alpha,beta,el);
    Cl = AeroCl(alpha,beta,el);

    delta_Cx_lef, delta_Cy_lef, delta_Cz_lef,delta_Cm_lef, delta_Cn_lef, delta_Cl_lef  = Delta_lef(alpha,beta); 

    # Aerodynamic Damping
    Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp  = AeroDamping(alpha);
    delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef, delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef, delta_Cnp_lef  = AeroDamping_lef(alpha);

    # Rudder Influence
    delta_Cy_r30, delta_Cn_r30, delta_Cl_r30 = AeroRudderInfluence(alpha,beta);

    # Aileron Influence
    delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef, delta_Cl_a20, delta_Cl_a20_lef = AeroAileronInfluence(alpha,beta);

    # Other Coefficients
    delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el = AeroOtherCoefficients(alpha,el);
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
  cbar = 11.32;       # mean aero chord, ft 
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
  Cx = AeroCx(alpha,beta,el); 
  Cy = AeroCy(alpha,beta);  
  Cz = AeroCz(alpha,beta,el);  
  Cm = AeroCm(alpha,beta,el);   
  Cn = AeroCn(alpha,beta,el);
  Cl = AeroCl(alpha,beta,el);

  delta_Cx_lef, delta_Cy_lef, delta_Cz_lef,delta_Cm_lef, delta_Cn_lef, delta_Cl_lef  = Delta_lef(alpha,beta); 

  # Aerodynamic Damping
  Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp  = AeroDamping(alpha);
  delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef, delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef, delta_Cnp_lef  = AeroDamping_lef(alpha);

  # Rudder Influence
  delta_Cy_r30, delta_Cn_r30, delta_Cl_r30 = AeroRudderInfluence(alpha,beta);

  # Aileron Influence
  delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef, delta_Cl_a20, delta_Cl_a20_lef = AeroAileronInfluence(alpha,beta);

  # Other Coefficients
  delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el = AeroOtherCoefficients(alpha,el);
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
  
