within TRANSFORM.Fluid.ClosureRelations.HeatTransfer.Models.DistributedPipe_1D_MultiTransferSurface;
partial model PartialTwoPhase
  extends PartialHeatTransfer_setQ_flows(
      redeclare replaceable package Medium =
        Modelica.Media.Water.StandardWater
      constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium);
  TRANSFORM.Media.BaseProperties2Phase[nHT] mediaProps(redeclare package Medium =
        Medium, state=states) "Bulk fluid properties"
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
//   TRANSFORM.Media.BaseProperties2Phase[nHT,nSurfaces] mediums_film(redeclare
//       package
//       Medium =
//         Medium, state=states_film) "Film fluid properties"
//     annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));
equation
  m_flows =vs .* mediaProps.d .* crossAreas;
  Res =mediaProps.d .* dimensions .* abs(vs) ./ mediaProps.mu;
  Prs =Medium.prandtlNumber(states);
//   for i in 1:nHT loop
//     for j in 1:nSurfaces loop
//       vs_film[i, j] = vs[i]*mediums[i].d/mediums_film[i,j].d;
//       Res_film[i, j] = mediums_film[i,j].d* dimensions[i]* abs(vs_film[i,j])/ mediums_film[i,j].mu;
//       Prs_film[i, j] = Medium.prandtlNumber(states_film[i,j]);
//     end for;
//   end for;
  annotation (defaultComponentName="heatTransfer",Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end PartialTwoPhase;
