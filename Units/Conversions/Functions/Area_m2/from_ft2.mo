within TRANSFORM.Units.Conversions.Functions.Area_m2;
function from_ft2 "Area: [ft2] -> [m2]"
  extends BaseClasses.from;
algorithm
  y := u*0.0254^2*12^2;
end from_ft2;
