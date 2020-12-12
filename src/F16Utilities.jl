# Utility files for the nonlinear F16 model.
#
# Created by:
# Raktim Bhattacharya,
# Aerospace Engineering, Texas A&M University
# isrlab.github.io
# March 22, 2020.
# --------------------------------------------------------------

using GridInterpolations, HDF5
 
struct GridData
    data::Vector;
    grid::RectangleGrid;
    nvar::Integer;
    minVal::Vector;
    maxVal::Vector; 
end

function CreateGridData(dataIn::Vector, inGrid::Vector...)::GridData
    data = dataIn;
    grid = RectangleGrid(inGrid...);
    nvar = length(inGrid);
    minVal = zeros(nvar);
    maxVal = zeros(nvar);
    for i=1:nvar
        minVal[i] = minimum(inGrid[i]);
        maxVal[i] = maximum(inGrid[i]);
    end
    return GridData(data,grid,nvar,minVal,maxVal);
end


# Interpolate with checks.
function evaluate(G::GridData,X::Vector; checkBounds::Bool=true)::Real
    if length(X) != G.nvar
        error("evaluate() --> Incompatible dimensions. Expecting vector of size $(G.nvar)");
    end

    if checkBounds
        outsideBound = false;
        for i=1:G.nvar
            if (X[i] > G.maxVal[i] || X[i] < G.minVal[i])
                    outsideBound = true;
                    println("evaluate() --> Extrapolating at index $i");
                    println("min: $(G.minVal[i]), val:$(X[i]), max: $(G.maxVal[i])");
            end
        end

        if outsideBound
            error("evaluate() --> Variable outside data range.");
        end
    end
    return interpolate(G.grid,G.data,X);
end


function readF16AeroData(fname::String)::Dict{String,GridData}
    # fid = h5open("F16AeroData.h5", "r");
    fid = h5open(fname, "r");

    # Read independent variables
    alpha1 = read(fid["alpha1"])[:];
    alpha2 = read(fid["alpha2"])[:];
    beta1 = read(fid["beta1"])[:];
    dh1 = read(fid["dh1"])[:];
    dh2 = read(fid["dh2"])[:];

    # Read Aero Tables and construct GridData data structures.
    F16AeroData = Dict{String,GridData}();

    # Add all the aero tables to the dictionary ... need to code this up.
    F16AeroData["_Cx"] = CreateGridData(read(fid["_Cx"])[:],alpha1,beta1,dh1);
    F16AeroData["_Cy"] = CreateGridData(read(fid["_Cy"])[:],alpha1,beta1);
    F16AeroData["_Cz"] = CreateGridData(read(fid["_Cz"])[:],alpha1,beta1,dh1);

    F16AeroData["_Cl"] = CreateGridData(read(fid["_Cl"])[:],alpha1,beta1,dh2);
    F16AeroData["_Cm"] = CreateGridData(read(fid["_Cm"])[:],alpha1,beta1,dh1);
    F16AeroData["_Cn"] = CreateGridData(read(fid["_Cn"])[:],alpha1,beta1,dh2);

    # Leading Edge Influence
    # ======================
    F16AeroData["_Cx_lef"] = CreateGridData(read(fid["_Cx_lef"])[:],alpha2,beta1);
    F16AeroData["_Cy_lef"] = CreateGridData(read(fid["_Cy_lef"])[:],alpha2,beta1);
    F16AeroData["_Cz_lef"] = CreateGridData(read(fid["_Cz_lef"])[:],alpha2,beta1);

    F16AeroData["_Cl_lef"] = CreateGridData(read(fid["_Cl_lef"])[:],alpha2,beta1);
    F16AeroData["_Cm_lef"] = CreateGridData(read(fid["_Cm_lef"])[:],alpha2,beta1);
    F16AeroData["_Cn_lef"] = CreateGridData(read(fid["_Cn_lef"])[:],alpha2,beta1);

    # Stability Derivatives
    # =====================
    F16AeroData["_Cxq"] = CreateGridData(read(fid["_Cxq"])[:],alpha1);
    F16AeroData["_Cyp"] = CreateGridData(read(fid["_Cyp"])[:],alpha1);
    F16AeroData["_Czq"] = CreateGridData(read(fid["_Czq"])[:],alpha1);
    F16AeroData["_Cmq"] = CreateGridData(read(fid["_Cmq"])[:],alpha1);

    F16AeroData["_Cyr"] = CreateGridData(read(fid["_Cyr"])[:],alpha1);
    F16AeroData["_Cnr"] = CreateGridData(read(fid["_Cnr"])[:],alpha1);

    F16AeroData["_Cnp"] = CreateGridData(read(fid["_Cnp"])[:],alpha1);
    F16AeroData["_Clp"] = CreateGridData(read(fid["_Clp"])[:],alpha1);
    F16AeroData["_Clr"] = CreateGridData(read(fid["_Clr"])[:],alpha1);


    F16AeroData["_deltaCxq_lef"] = CreateGridData(read(fid["_deltaCxq_lef"])[:],alpha2);
    F16AeroData["_deltaCyr_lef"] = CreateGridData(read(fid["_deltaCyr_lef"])[:],alpha2);
    F16AeroData["_deltaCyp_lef"] = CreateGridData(read(fid["_deltaCyp_lef"])[:],alpha2);

    F16AeroData["_deltaCzq_lef"] = CreateGridData(read(fid["_deltaCzq_lef"])[:],alpha2);
    F16AeroData["_deltaClr_lef"] = CreateGridData(read(fid["_deltaClr_lef"])[:],alpha2);
    F16AeroData["_deltaClp_lef"] = CreateGridData(read(fid["_deltaClp_lef"])[:],alpha2);

    F16AeroData["_deltaCmq_lef"] = CreateGridData(read(fid["_deltaCmq_lef"])[:],alpha2);
    F16AeroData["_deltaCnr_lef"] = CreateGridData(read(fid["_deltaCnr_lef"])[:],alpha2);
    F16AeroData["_deltaCnp_lef"] = CreateGridData(read(fid["_deltaCnp_lef"])[:],alpha2);

    F16AeroData["_Cy_r30"] = CreateGridData(read(fid["_Cy_r30"])[:],alpha1,beta1);
    F16AeroData["_Cn_r30"] = CreateGridData(read(fid["_Cn_r30"])[:],alpha1,beta1);
    F16AeroData["_Cl_r30"] = CreateGridData(read(fid["_Cl_r30"])[:],alpha1,beta1);

    F16AeroData["_Cy_a20"] = CreateGridData(read(fid["_Cy_a20"])[:],alpha1,beta1);
    F16AeroData["_Cy_a20_lef"] = CreateGridData(read(fid["_Cy_a20_lef"])[:],alpha2,beta1);

    F16AeroData["_Cn_a20"] = CreateGridData(read(fid["_Cn_a20"])[:],alpha1,beta1);
    F16AeroData["_Cn_a20_lef"] = CreateGridData(read(fid["_Cn_a20_lef"])[:],alpha2,beta1);

    F16AeroData["_Cl_a20"] = CreateGridData(read(fid["_Cl_a20"])[:],alpha1,beta1);
    F16AeroData["_Cl_a20_lef"] = CreateGridData(read(fid["_Cl_a20_lef"])[:],alpha2,beta1);

    F16AeroData["_deltaCnbeta"] = CreateGridData(read(fid["_deltaCnbeta"])[:],alpha1);
    F16AeroData["_deltaClbeta"] = CreateGridData(read(fid["_deltaClbeta"])[:],alpha1);
    F16AeroData["_deltaCm"] = CreateGridData(read(fid["_deltaCm"])[:],alpha1);

    F16AeroData["_eta_el"] = CreateGridData(read(fid["_eta_el"])[:],dh1);

    close(fid);

    return F16AeroData;
end
