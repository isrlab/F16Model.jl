module F16Model
using FileIO, GridInterpolations, ForwardDiff

global F16AeroData = load("data/F16AeroData.jld2"); # Loads the Aero Tables for F16.


#include all the files
include("NonlinearF16Model.jl");

end # module
