within TRANSFORM.Fluid.Machines.BaseClasses.TurbineCharacteristics.Efficiency;
model Constant "Constant efficiency"
  extends PartialEfficiencyChar;

  parameter SI.Efficiency eta_constant=1.0 "Constant efficiency"
    annotation (Dialog);

equation
  eta =  eta_constant;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Constant;
