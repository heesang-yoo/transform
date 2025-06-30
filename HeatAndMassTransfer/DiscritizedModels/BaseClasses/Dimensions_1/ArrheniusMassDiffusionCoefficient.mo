within TRANSFORM.HeatAndMassTransfer.DiscritizedModels.BaseClasses.Dimensions_1;
model ArrheniusMassDiffusionCoefficient
  import TRANSFORM.Math.fillArray_1D;
  extends PartialMassDiffusionCoefficient;
  input SI.Temperature T[nVs[1]] = Material.temperature(states) "Temperature" annotation(Dialog(group="Inputs"));
  input SI.MolarEnergy Ea[nVs[1],nC] "Activation energy" annotation(Dialog(group="Inputs"));
  input Real A[nVs[1],nC] "Pre-exponential factor" annotation(Dialog(group="Inputs"));
  input Real beta[nVs[1],nC] = fill(1.0,nVs[1],nC) "Correction factor" annotation(Dialog(group="Inputs"));
  input SI.MolarHeatCapacity R = Modelica.Constants.R "Universal gas constant" annotation(Dialog(group="Inputs"));
equation
  for ic in 1:nC loop
  for i in 1:nVs[1] loop
    D_abs[i,ic] = A[i,ic]* exp(-(Ea[i,ic]/ (R* T[i]))^ beta[i,ic]);
  end for;
  end for;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ArrheniusMassDiffusionCoefficient;
