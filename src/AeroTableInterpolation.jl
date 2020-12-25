using HDF5, Interpolations

fname = joinpath(dirname(pathof(F16Model)), "data/F16AeroData.h5");
fid = h5open(fname, "r");

function createAeroFunction(data,indepVars...)
    len = length.(indepVars);
    y = reshape(data,len);
    itp = interpolate(indepVars,y,Gridded(Linear()));
    return itp;
end

# Read independent variables
alpha1 = read(fid["alpha1"])[:];
alpha2 = read(fid["alpha2"])[:];
beta1 = read(fid["beta1"])[:];
dh1 = read(fid["dh1"])[:];
dh2 = read(fid["dh2"])[:];

# Create aero functions
_Cx = createAeroFunction(read(fid["_Cx"])[:],alpha1,beta1,dh1);
_Cy = createAeroFunction(read(fid["_Cy"])[:],alpha1,beta1);
_Cz = createAeroFunction(read(fid["_Cz"])[:],alpha1,beta1,dh1);

_Cl = createAeroFunction(read(fid["_Cl"])[:],alpha1,beta1,dh2);
_Cm = createAeroFunction(read(fid["_Cm"])[:],alpha1,beta1,dh1);
_Cn = createAeroFunction(read(fid["_Cn"])[:],alpha1,beta1,dh2);

# Leading Edge Influence
# ======================
_Cx_lef = createAeroFunction(read(fid["_Cx_lef"])[:],alpha2,beta1);
_Cy_lef = createAeroFunction(read(fid["_Cy_lef"])[:],alpha2,beta1);
_Cz_lef = createAeroFunction(read(fid["_Cz_lef"])[:],alpha2,beta1);

_Cl_lef = createAeroFunction(read(fid["_Cl_lef"])[:],alpha2,beta1);
_Cm_lef = createAeroFunction(read(fid["_Cm_lef"])[:],alpha2,beta1);
_Cn_lef = createAeroFunction(read(fid["_Cn_lef"])[:],alpha2,beta1);

# Stability Derivatives
# =====================
_Cxq = createAeroFunction(read(fid["_Cxq"])[:],alpha1);
_Cyp = createAeroFunction(read(fid["_Cyp"])[:],alpha1);
_Czq = createAeroFunction(read(fid["_Czq"])[:],alpha1);
_Cmq = createAeroFunction(read(fid["_Cmq"])[:],alpha1);

_Cyr = createAeroFunction(read(fid["_Cyr"])[:],alpha1);
_Cnr = createAeroFunction(read(fid["_Cnr"])[:],alpha1);

_Cnp = createAeroFunction(read(fid["_Cnp"])[:],alpha1);
_Clp = createAeroFunction(read(fid["_Clp"])[:],alpha1);
_Clr = createAeroFunction(read(fid["_Clr"])[:],alpha1);

_deltaCxq_lef = createAeroFunction(read(fid["_deltaCxq_lef"])[:],alpha2);
_deltaCyr_lef = createAeroFunction(read(fid["_deltaCyr_lef"])[:],alpha2);
_deltaCyp_lef = createAeroFunction(read(fid["_deltaCyp_lef"])[:],alpha2);

_deltaCzq_lef = createAeroFunction(read(fid["_deltaCzq_lef"])[:],alpha2);
_deltaClr_lef = createAeroFunction(read(fid["_deltaClr_lef"])[:],alpha2);
_deltaClp_lef = createAeroFunction(read(fid["_deltaClp_lef"])[:],alpha2);

_deltaCmq_lef = createAeroFunction(read(fid["_deltaCmq_lef"])[:],alpha2);
_deltaCnr_lef = createAeroFunction(read(fid["_deltaCnr_lef"])[:],alpha2);
_deltaCnp_lef = createAeroFunction(read(fid["_deltaCnp_lef"])[:],alpha2);

_Cy_r30 = createAeroFunction(read(fid["_Cy_r30"])[:],alpha1,beta1);
_Cn_r30 = createAeroFunction(read(fid["_Cn_r30"])[:],alpha1,beta1);
_Cl_r30 = createAeroFunction(read(fid["_Cl_r30"])[:],alpha1,beta1);

_Cy_a20 = createAeroFunction(read(fid["_Cy_a20"])[:],alpha1,beta1);
_Cy_a20_lef = createAeroFunction(read(fid["_Cy_a20_lef"])[:],alpha2,beta1);

_Cn_a20 = createAeroFunction(read(fid["_Cn_a20"])[:],alpha1,beta1);
_Cn_a20_lef = createAeroFunction(read(fid["_Cn_a20_lef"])[:],alpha2,beta1);

_Cl_a20 = createAeroFunction(read(fid["_Cl_a20"])[:],alpha1,beta1);
_Cl_a20_lef = createAeroFunction(read(fid["_Cl_a20_lef"])[:],alpha2,beta1);

_deltaCnbeta = createAeroFunction(read(fid["_deltaCnbeta"])[:],alpha1);
_deltaClbeta = createAeroFunction(read(fid["_deltaClbeta"])[:],alpha1);
_deltaCm = createAeroFunction(read(fid["_deltaCm"])[:],alpha1);

_eta_el = createAeroFunction(read(fid["_eta_el"])[:],dh1);

close(fid);


# ------------- Other derived aero functions ----------------

function Delta_lef(alpha::Real,beta::Real)::Tuple{Real,Real,Real,Real,Real,Real}
    delta_Cx_lef = _Cx_lef(alpha,beta) - _Cx(alpha,beta,0.0);
    delta_Cy_lef = _Cy_lef(alpha,beta) - _Cy(alpha,beta);
    delta_Cz_lef = _Cz_lef(alpha,beta) - _Cz(alpha,beta,0.0);

    delta_Cl_lef = _Cl_lef(alpha,beta) - _Cl(alpha,beta,0.0);
    delta_Cm_lef = _Cm_lef(alpha,beta) - _Cm(alpha,beta,0.0);
    delta_Cn_lef = _Cn_lef(alpha,beta) - _Cn(alpha,beta,0.0);

    return (delta_Cx_lef, delta_Cy_lef, delta_Cz_lef,delta_Cm_lef, delta_Cn_lef, delta_Cl_lef);
end

# Damping Terms
function _Damping(alpha::Real)::Tuple{Real,Real,Real,Real,Real,Real,Real,Real,Real}
    Cxq = _Cxq(alpha);
    Cyr = _Cyr(alpha);
    Cyp = _Cyp(alpha);
    Czq = _Czq(alpha);
    Clr = _Clr(alpha);
    Clp = _Clp(alpha);
    Cmq = _Cmq(alpha);
    Cnr = _Cnr(alpha);
    Cnp = _Cnp(alpha);
    
    return (Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp);
end

function _Damping_lef(alpha::Real)::Tuple{Real,Real,Real,Real,Real,Real,Real,Real,Real}
    delta_Cxq_lef = _deltaCxq_lef(alpha);
    delta_Cyr_lef = _deltaCyr_lef(alpha);
    delta_Cyp_lef = _deltaCyp_lef(alpha);
    delta_Czq_lef = _deltaCzq_lef(alpha);
    delta_Clr_lef = _deltaClr_lef(alpha);
    delta_Clp_lef = _deltaClp_lef(alpha);
    delta_Cmq_lef = _deltaCmq_lef(alpha);
    delta_Cnr_lef = _deltaCnr_lef(alpha);
    delta_Cnp_lef = _deltaCnp_lef(alpha);
    
    return (delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef, delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef, delta_Cnp_lef);
end


# Rudder Influence
function _RudderInfluence(alpha::Real,beta::Real)::Tuple{Real,Real,Real}
    delta_Cy_r30 = _Cy_r30(alpha,beta) - _Cy(alpha,beta);
    delta_Cn_r30 = _Cn_r30(alpha,beta) - _Cn(alpha,beta,0.0);
    delta_Cl_r30 = _Cl_r30(alpha,beta) - _Cl(alpha,beta,0.0);

    return (delta_Cy_r30, delta_Cn_r30, delta_Cl_r30);
end

 function _AileronInfluence(alpha::Real,beta::Real)::Tuple{Real,Real,Real,Real,Real,Real}
    delta_Cy_a20 = _Cy_a20(alpha,beta) - _Cy(alpha,beta);
    delta_Cy_a20_lef = _Cy_a20_lef(alpha,beta) - _Cy_lef(alpha,beta) - delta_Cy_a20;

    delta_Cn_a20 = _Cn_a20(alpha,beta) - _Cn(alpha,beta,0.0);
    delta_Cn_a20_lef = _Cn_a20_lef(alpha,beta) - _Cn_lef(alpha,beta) - delta_Cn_a20;

    delta_Cl_a20 = _Cl_a20(alpha,beta) - _Cl(alpha,beta,0.0);
    delta_Cl_a20_lef = _Cl_a20_lef(alpha,beta) - _Cl_lef(alpha,beta) - delta_Cl_a20;

    return (delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef, delta_Cl_a20, delta_Cl_a20_lef);
 end

 function _OtherCoefficients(alpha,el)::Tuple{Real,Real,Real,Real}
    delta_Cnbeta = _deltaCnbeta(alpha);
    delta_Clbeta = _deltaClbeta(alpha);
    delta_Cm = _deltaCm(alpha);
    eta_el = _eta_el(el);

    return (delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el);
 end
