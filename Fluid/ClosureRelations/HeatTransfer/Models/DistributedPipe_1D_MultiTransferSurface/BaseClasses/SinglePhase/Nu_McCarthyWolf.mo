within TRANSFORM.Fluid.ClosureRelations.HeatTransfer.Models.DistributedPipe_1D_MultiTransferSurface.BaseClasses.SinglePhase;
model Nu_McCarthyWolf "McCarthy-Wolf"
  extends PartialHeatTransferCorrelation;

  input Real A =  0.025 "Multiplication value" annotation(Dialog(group="Inputs"));
  input Real alpha =  0.8   "Exponent to Reynolds number" annotation(Dialog(group="Inputs"));
  input Real beta =  0.4   "Exponent to Prandtl number" annotation(Dialog(group="Inputs"));
  input Real gamma = -0.55  "Exponent for temperature ratio" annotation(Dialog(group="Inputs"));

  SI.Temperature T_wall = Medium.temperature(state_wall);

equation
  Nu = A*Re^alpha*Pr^beta*(T_wall/mediaProps.T)^gamma;

  annotation (defaultComponentName="heatTransfer",Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Nu_McCarthyWolf;
