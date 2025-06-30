within TRANSFORM.HeatExchangers.ClosureRelations.Models.EffectivenessNTU_Relations;
partial model PartialMethod
  extends PartialEffectivenessNTU;
  Units.NonDim method "Method for C_r > 0";
  parameter Boolean epsilonMethod=true
    "=true to use formulation of epsilon=f(NTU) else NTU=f(epsilon)"
    annotation (Dialog(tab="Advanced"));
equation
  if epsilonMethod then
    epsilon = Math.spliceTanh(
      method,
      1 - exp(-NTU),
      C_r,
      0.01);
  else
    NTU = Math.spliceTanh(
      method,
      -log(1 - epsilon),
      C_r,
      0.01);
  end if;
  annotation (
    defaultComponentName="effectivenessNTU",
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)));
end PartialMethod;
