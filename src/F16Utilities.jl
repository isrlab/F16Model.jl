# Utility files for the nonlinear F16 model.
#
# Created by:
# Raktim Bhattacharya,
# Aerospace Engineering, Texas A&M University
# isrlab.github.io
# March 22, 2020.
# --------------------------------------------------------------

using GridInterpolations, LinearAlgebra
 
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
            error("evaluate() --> Variable outside _ range.");
        end
    end
    return interpolate(G.grid,G.data,X);
end