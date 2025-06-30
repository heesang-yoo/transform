within TRANSFORM.Fluid.ClosureRelations.MassTransfer.Models.DistributedPipe_TraceMass_1D_MultiTransferSurface;
model AlphasM "Specify Mass Transfer Coefficient (alphaM)"
  import TRANSFORM.Math.fillArray_2D;
  extends PartialSinglePhase;
  input Units.CoefficientOfMassTransfer alphaM0[nC]=fill(0, nC)
    "Coefficient of mass transfer" annotation (Dialog(group="Inputs"));
  input Units.CoefficientOfMassTransfer alphasM0[nMT,nSurfaces,nC]=fillArray_2D(
      alphaM0, nMT,nSurfaces) "if non-uniform then set"
    annotation (Dialog(group="Inputs"));
equation
  for i in 1:nMT loop
    for j in 1:nSurfaces loop
      alphasM[i, j, :] = alphasM0[i, j, :];
      Shs[i, j, :] = alphasM[i, j, :] .* dimensions[i] ./ diffusionCoeff[i].D_abs;
    end for;
  end for;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end AlphasM;
