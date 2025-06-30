within TRANSFORM.Fluid.ClosureRelations.PressureLoss.Models.DistributedPipe_1D;
model SinglePhase_ReynoldsRelation_NumStable
  "Single Phase | Reynolds Relation | Numerically Stable Method"
  extends PartialSinglePhase;
  TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_IN_con[
    nFM] IN_con(
    diameter_a=dimensions[1:nFM],
    diameter_b=dimensions[2:nFM + 1],
    crossArea_a=crossAreas[1:nFM],
    crossArea_b=crossAreas[2:nFM + 1],
    each A=A,
    each B=B,
    each C=C)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_IN_var[
    nFM] IN_var(
    rho_a=ds_a,
    rho_b=ds_b,
    mu_a=mus_a,
    mu_b=mus_b)
    annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
  //Reynolds Relation parameters
  parameter Real A;
  parameter Real B;
  parameter Real C;
  // Parameters
  parameter SI.AbsolutePressure dp_nominal(start=1, fixed=false)
    "Nominal pressure loss (only for nominal models)" annotation (Dialog(tab="Advanced", group="Nominal Conditions", enable=false));
  parameter SI.MassFlowRate m_flow_nominal=1e2*m_flow_small
    "Nominal mass flow rate"
    annotation (Dialog(tab="Advanced", group="Nominal Conditions"));
  parameter Boolean use_d_nominal=false
    "= true, if d_nominal is used, otherwise computed from medium" annotation (
      Dialog(tab="Advanced", group="Nominal Conditions"), Evaluate=true);
  parameter SI.Density d_nominal=Medium.density_phX(
      Medium.p_default,
      Medium.h_default,
      Medium.X_default)
    "Nominal density (e.g., rho_liquidWater = 995, rho_air = 1.2)" annotation (
      Dialog(
      tab="Advanced",
      group="Nominal Conditions",
      enable=use_d_nominal));
  parameter Boolean use_mu_nominal=false
    "= true, if mu_nominal is used, otherwise computed from medium" annotation (
     Dialog(tab="Advanced", group="Nominal Conditions"), Evaluate=true);
  parameter SI.DynamicViscosity mu_nominal=Medium.dynamicViscosity(
      Medium.setState_phX(
      Medium.p_default,
      Medium.h_default,
      Medium.X_default))
    "Nominal dynamic viscosity (e.g., mu_liquidWater = 1e-3, mu_air = 1.8e-5)"
    annotation (Dialog(
      tab="Advanced",
      group="Nominal Conditions",
      enable=use_mu_nominal));
  parameter Boolean continuousFlowReversal=if use_d_nominal and use_mu_nominal
       then true else false
    "= true if the pressure loss is continuous around zero flow" annotation (
      Dialog(tab="Advanced", group="Nominal Conditions"), Evaluate=true);
  SI.DynamicViscosity mus[nFM];
  SI.DynamicViscosity mus_a[nFM]=if use_mu_nominal then fill(mu_nominal, nFM)
       else mediaProps[1:nFM].mu;
  SI.DynamicViscosity mus_b[nFM]=if use_mu_nominal then fill(mu_nominal, nFM)
       else mediaProps[2:nFM + 1].mu;
  SI.Density ds[nFM];
  SI.Density ds_a[nFM]=if use_d_nominal then fill(d_nominal, nFM) else
      mediaProps[1:nFM].d;
  SI.Density ds_b[nFM]=if use_d_nominal then fill(d_nominal, nFM) else
      mediaProps[2:nFM + 1].d;
protected
  TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_IN_var
    IN_var_nominal(
    rho_a=d_nominal,
    rho_b=d_nominal,
    mu_a=mu_nominal,
    mu_b=mu_nominal);
  SI.AbsolutePressure dp_fric_nominal=sum(
      TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_DP(
      IN_con,
      IN_var_nominal,
      m_flow_nominal,
      m_flow_small)) "pressure loss for nominal conditions";
initial equation
  dp_nominal = dp_fric_nominal + g_n*sum(dheights)*d_nominal;
equation
  for i in 1:nFM loop
    ds[i] = TRANSFORM.Math.spliceTanh(
      ds_a[i],
      ds_b[i],
      m_flows[i],
      m_flow_small);
    mus[i] = TRANSFORM.Math.spliceTanh(
      mus_a[i],
      mus_b[i],
      m_flows[i],
      m_flow_small);
  end for;
  if continuousFlowReversal then
    // simple regularization
    if from_dp then
      // and not dp_is_zero then
      m_flows =homotopy(actual=
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_MFLOW(
        IN_con,
        IN_var,
        dps_fg - {g_n*dheights[i]*ds_a[i] for i in 1:nFM},
        dp_small/(nFM)), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
        dheights*d_nominal));
    else
      dps_fg =homotopy(actual=
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_DP(
        IN_con,
        IN_var,
        m_flows,
        m_flow_small) + {g_n*dheights[i]*ds_a[i] for i in 1:nFM}, simplified=
        dp_nominal/m_flow_nominal*m_flows + g_n*dheights*d_nominal);
    end if;
  else
    // regularization for discontinuous flow reversal and static head
    if from_dp then
      // and not dp_is_zero then
      m_flows =homotopy(actual=
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_MFLOW_staticHead(
        IN_con,
        IN_var,
        dps_fg,
        dp_small/(nFM),
        g_n*dheights), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
        dheights*d_nominal));
    else
      dps_fg =homotopy(actual=
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_ReynoldsRelation.dp_DP_staticHead(
        IN_con,
        IN_var,
        m_flows,
        m_flow_small,
        g_n*dheights), simplified=dp_nominal/m_flow_nominal*m_flows + g_n*
        dheights*d_nominal);
    end if;
  end if;
    annotation (Dialog(tab="Advanced", group="Nominal Conditions",enable=false),
              Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end SinglePhase_ReynoldsRelation_NumStable;
