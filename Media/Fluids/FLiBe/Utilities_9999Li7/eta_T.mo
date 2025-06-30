within TRANSFORM.Media.Fluids.FLiBe.Utilities_9999Li7;
function eta_T
  import TRANSFORM.Units.Conversions.Functions.Temperature_K.to_degR;
  import from_lb_hrfeet =
    TRANSFORM.Units.Conversions.Functions.DynamicViscosity_kg_ms.from_lb_hrft;
  input SI.Temperature T;
  output SI.DynamicViscosity eta;
algorithm
  eta:=from_lb_hrfeet(0.2806*exp(6759/to_degR(T)));
end eta_T;
