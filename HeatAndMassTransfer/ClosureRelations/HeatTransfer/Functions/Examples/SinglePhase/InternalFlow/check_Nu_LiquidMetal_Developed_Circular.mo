within TRANSFORM.HeatAndMassTransfer.ClosureRelations.HeatTransfer.Functions.Examples.SinglePhase.InternalFlow;
model check_Nu_LiquidMetal_Developed_Circular
  extends TRANSFORM.Icons.Example;
  parameter SI.PrandtlNumber Pr = 1.5 "Prandtl Number";
  parameter Boolean constantTwall=true
    "= true for constant wall temperature correlation else constant heat flux (laminar conditions only Re ~< 2300)";
    Real y[2];
  Modelica.Blocks.Sources.Trapezoid Re(
    rising=1,
    width=0,
    falling=1,
    period=2,
    nperiod=1,
    amplitude=4000)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  Utilities.ErrorAnalysis.UnitTests unitTests(n=2, x=y)
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
equation
 y[1] = TRANSFORM.HeatAndMassTransfer.ClosureRelations.HeatTransfer.Functions.SinglePhase.InternalFlow.Nu_LiquidMetal_Developed_Circular(Re.y, Pr, constantTwall);
 y[2] = TRANSFORM.HeatAndMassTransfer.ClosureRelations.HeatTransfer.Functions.SinglePhase.InternalFlow.Nu_LiquidMetal_Developed_Circular(Re.y, Pr, not constantTwall);
  annotation (experiment(StopTime=2, __Dymola_NumberOfIntervals=160));
end check_Nu_LiquidMetal_Developed_Circular;
