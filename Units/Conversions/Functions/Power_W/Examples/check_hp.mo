within TRANSFORM.Units.Conversions.Functions.Power_W.Examples;
model check_hp
  extends TRANSFORM.Icons.Example;
  parameter SI.Length u=1;
  final parameter Real x_reference[unitTests.n]={0.00134102,1/0.00134102};
  Real x[unitTests.n];
  Utilities.ErrorAnalysis.UnitTests unitTests(n=2,x=x, x_reference=x_reference)
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
  x[1] = to_hp(u);
  x[2] = from_hp(u);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end check_hp;
