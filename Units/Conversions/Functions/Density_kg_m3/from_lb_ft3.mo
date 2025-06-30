within TRANSFORM.Units.Conversions.Functions.Density_kg_m3;
function from_lb_ft3 "Density: [lb/ft^3] -> [kg/m^3]"
  extends BaseClasses.from;
algorithm
  y := u*0.453592/(12^3*0.0254^3);
end from_lb_ft3;
