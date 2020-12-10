
using FileIO

include("F16Utilities.jl")
global F16AeroData = load("F16AeroData.jld2"); # Loads the Aero Tables for F16.

# ------------- Code for Table Interpolation ----------------

# Forces
function AeroCx(alpha::Real,beta::Real,el::Real)::Real
    return evaluate(F16AeroData["_Cx"],[alpha, beta, el]);
end

function AeroCy(alpha::Real,beta::Real)::Real
    return evaluate(F16AeroData["_Cy"],[alpha, beta]);
end

function AeroCz(alpha::Real,beta::Real,el::Real)::Real
    return evaluate(F16AeroData["_Cz"],[alpha, beta, el]);
end

# Moments
function AeroCl(alpha::Real,beta::Real,el::Real)::Real
    return evaluate(F16AeroData["_Cl"],[alpha, beta, el]);
end

function AeroCm(alpha::Real,beta::Real,el::Real)::Real
    return evaluate(F16AeroData["_Cm"],[alpha, beta, el]);
end

function AeroCn(alpha::Real,beta::Real,el::Real)::Real
    return evaluate(F16AeroData["_Cn"],[alpha, beta, el]);
end

# Leading Edge Influence
function AeroCx_lef(alpha::Real,beta::Real)
    return evaluate(F16AeroData["_Cx_lef"],[alpha, beta]);
end

function AeroCy_lef(alpha::Real,beta::Real)
    return evaluate(F16AeroData["_Cy_lef"],[alpha, beta]);
end

function AeroCz_lef(alpha::Real,beta::Real)
    return evaluate(F16AeroData["_Cz_lef"],[alpha, beta]);
end

function AeroCl_lef(alpha::Real,beta::Real)
    return evaluate(F16AeroData["_Cl_lef"],[alpha, beta]);
end

function AeroCm_lef(alpha::Real,beta::Real)
    return evaluate(F16AeroData["_Cm_lef"],[alpha, beta]);
end

function AeroCn_lef(alpha::Real,beta::Real)
    return evaluate(F16AeroData["_Cn_lef"],[alpha, beta]);
end

function Delta_lef(alpha::Real,beta::Real)::Tuple{Real,Real,Real,Real,Real,Real}
    delta_Cx_lef = AeroCx_lef(alpha,beta) - AeroCx(alpha,beta,0.0);
    delta_Cy_lef = AeroCy_lef(alpha,beta) - AeroCy(alpha,beta);
    delta_Cz_lef = AeroCz_lef(alpha,beta) - AeroCz(alpha,beta,0.0);

    delta_Cl_lef = AeroCl_lef(alpha,beta) - AeroCl(alpha,beta,0.0);
    delta_Cm_lef = AeroCm_lef(alpha,beta) - AeroCm(alpha,beta,0.0);
    delta_Cn_lef = AeroCn_lef(alpha,beta) - AeroCn(alpha,beta,0.0);

    return (delta_Cx_lef, delta_Cy_lef, delta_Cz_lef,delta_Cm_lef, delta_Cn_lef, delta_Cl_lef);
end

# Damping Terms
function AeroDamping(alpha::Real)::Tuple{Real,Real,Real,Real,Real,Real,Real,Real,Real}
    Cxq = evaluate(F16AeroData["_Cxq"],[alpha]);
    Cyr = evaluate(F16AeroData["_Cyr"],[alpha]);
    Cyp = evaluate(F16AeroData["_Cyp"],[alpha]);
    Czq = evaluate(F16AeroData["_Czq"],[alpha]);
    Clr = evaluate(F16AeroData["_Clr"],[alpha]);
    Clp = evaluate(F16AeroData["_Clp"],[alpha]);
    Cmq = evaluate(F16AeroData["_Cmq"],[alpha]);
    Cnr = evaluate(F16AeroData["_Cnr"],[alpha]);
    Cnp = evaluate(F16AeroData["_Cnp"],[alpha]);
    
    return (Cxq, Cyr, Cyp, Czq, Clr, Clp, Cmq, Cnr, Cnp);
end

function AeroDamping_lef(alpha::Real)::Tuple{Real,Real,Real,Real,Real,Real,Real,Real,Real}
    delta_Cxq_lef = evaluate(F16AeroData["_deltaCxq_lef"],[alpha]);
    delta_Cyr_lef = evaluate(F16AeroData["_deltaCyr_lef"],[alpha]);
    delta_Cyp_lef = evaluate(F16AeroData["_deltaCyp_lef"],[alpha]);
    delta_Czq_lef = evaluate(F16AeroData["_deltaCzq_lef"],[alpha]);
    delta_Clr_lef = evaluate(F16AeroData["_deltaClr_lef"],[alpha]);
    delta_Clp_lef = evaluate(F16AeroData["_deltaClp_lef"],[alpha]);
    delta_Cmq_lef = evaluate(F16AeroData["_deltaCmq_lef"],[alpha]);
    delta_Cnr_lef = evaluate(F16AeroData["_deltaCnr_lef"],[alpha]);
    delta_Cnp_lef = evaluate(F16AeroData["_deltaCnp_lef"],[alpha]);
    
    return (delta_Cxq_lef, delta_Cyr_lef, delta_Cyp_lef, delta_Czq_lef, delta_Clr_lef, delta_Clp_lef, delta_Cmq_lef, delta_Cnr_lef, delta_Cnp_lef);
end


# Rudder Influence
function AeroRudderInfluence(alpha::Real,beta::Real)::Tuple{Real,Real,Real}
    delta_Cy_r30 = evaluate(F16AeroData["_Cy_r30"],[alpha,beta]) - AeroCy(alpha,beta);
    delta_Cn_r30 = evaluate(F16AeroData["_Cn_r30"],[alpha,beta]) - AeroCn(alpha,beta,0.0);
    delta_Cl_r30 = evaluate(F16AeroData["_Cl_r30"],[alpha,beta]) - AeroCl(alpha,beta,0.0);

    return (delta_Cy_r30, delta_Cn_r30, delta_Cl_r30);
end

 function AeroAileronInfluence(alpha::Real,beta::Real)::Tuple{Real,Real,Real,Real,Real,Real}
    delta_Cy_a20 = evaluate(F16AeroData["_Cy_a20"],[alpha,beta]) - AeroCy(alpha,beta);
    delta_Cy_a20_lef = evaluate(F16AeroData["_Cy_a20_lef"],[alpha,beta]) - AeroCy_lef(alpha,beta) - delta_Cy_a20;

    delta_Cn_a20 = evaluate(F16AeroData["_Cn_a20"],[alpha,beta]) - AeroCn(alpha,beta,0.0);
    delta_Cn_a20_lef = evaluate(F16AeroData["_Cn_a20_lef"],[alpha,beta]) - AeroCn_lef(alpha,beta) - delta_Cn_a20;

    delta_Cl_a20 = evaluate(F16AeroData["_Cl_a20"],[alpha,beta]) - AeroCl(alpha,beta,0.0);
    delta_Cl_a20_lef = evaluate(F16AeroData["_Cl_a20_lef"],[alpha,beta]) - AeroCl_lef(alpha,beta) - delta_Cl_a20;

    return (delta_Cy_a20, delta_Cy_a20_lef, delta_Cn_a20, delta_Cn_a20_lef, delta_Cl_a20, delta_Cl_a20_lef);
 end

 function AeroOtherCoefficients(alpha,el)::Tuple{Real,Real,Real,Real}
    delta_Cnbeta = evaluate(F16AeroData["_deltaCnbeta"],[alpha]);
    delta_Clbeta = evaluate(F16AeroData["_deltaClbeta"],[alpha]);
    delta_Cm = evaluate(F16AeroData["_deltaCm"],[alpha]);
    eta_el = evaluate(F16AeroData["_eta_el"],[el]);

    return (delta_Cnbeta, delta_Clbeta, delta_Cm, eta_el);
 end
