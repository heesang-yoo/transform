within TRANSFORM.Media.Fluids.FLiBe.Utilities_12Th_05U;
function cp_T
  import from_btu_lbF =
    TRANSFORM.Units.Conversions.Functions.SpecificHeatAndEntropy_J_kgK.from_btu_lbdegF;
  input SI.Temperature T;
  output SI.SpecificHeatCapacity cp;
algorithm
  cp:=from_btu_lbF(0.32);
end cp_T;
