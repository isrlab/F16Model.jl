module F16Model

using GridInterpolations, ForwardDiff, HDF5

include("F16Utilities.jl");

fname = joinpath(dirname(pathof(F16Model)), "F16AeroData.h5");
global F16AeroData = readF16AeroData(fname);

#include all the files
include("NonlinearF16Model.jl");

end # module
