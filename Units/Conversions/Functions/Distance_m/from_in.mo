within TRANSFORM.Units.Conversions.Functions.Distance_m;
function from_in "Distance: [in] -> [m]"
  extends BaseClasses.from;
algorithm
  y := u*0.0254;
end from_in;
