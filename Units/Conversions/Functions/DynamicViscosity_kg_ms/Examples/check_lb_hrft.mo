within TRANSFORM.Units.Conversions.Functions.DynamicViscosity_kg_ms.Examples;
model check_lb_hrft
  extends TRANSFORM.Icons.Example;
  parameter SI.Length u=1;
  final parameter Real x_reference[unitTests.n]={2419.0883293091,1/2419.0883293091};
  Real x[unitTests.n];
  Utilities.ErrorAnalysis.UnitTests unitTests(n=2,x=x, x_reference=x_reference)
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
  x[1] = to_lb_hrft(u);
  x[2] = from_lb_hrft(u);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end check_lb_hrft;
