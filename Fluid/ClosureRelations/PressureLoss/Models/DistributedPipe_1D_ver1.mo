within TRANSFORM.Fluid.ClosureRelations.PressureLoss.Models;
package DistributedPipe_1D_ver1
  extends Icons.VariantsPackage;
  partial model PartialDistributedStaggeredFlow
    import Modelica.Fluid.Types.Dynamics;
    replaceable package Medium = Modelica.Media.Water.StandardWater
      constrainedby Modelica.Media.Interfaces.PartialMedium "Medium properties"
      annotation (choicesAllMatching=true, Dialog(tab="Internal Interface"));
    parameter Integer nFM(min=1) = 1 "Number of discrete flow models"
      annotation (Dialog(tab="Internal Interface"));
    input SI.Acceleration g_n=Modelica.Constants.g_n "Gravitational acceleration"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    replaceable input SI.Acceleration g_eff[nFM] = fill(Modelica.Constants.g_n, nFM)  "Effective Gravitational acceleration";
    // Initialization
    parameter Dynamics momentumDynamics = Dynamics.DynamicFreeInitial "Formulation of momentum balance"
      annotation (Dialog(tab="Internal Interface", group="Initialization"));
    parameter SI.PressureDifference[nFM] dps_start
      "Pressure changes {p[2]-p[1],...,p[n+1]-p[n]}"
      annotation (Dialog(tab="Internal Interface", group="Initialization"));
    parameter SI.MassFlowRate[nFM] m_flows_start "Mass flow rates"
      annotation (Dialog(tab="Internal Interface", group="Initialization"));
    // State parameters
    input Medium.ThermodynamicState states[nFM + 1] "State at volumes"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Velocity[nFM + 1] vs "Mean velocities of fluid flow"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Temperature[nFM + 1] Ts_wall
      "Mean wall temperatures of heat transfer surface"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    // Geometry
    input SI.Length dimensions[nFM + 1]
      "Characteristic dimensions (e.g. hydraulic diameter)"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Area crossAreas[nFM + 1] "Cross sectional area"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Length perimeters[nFM + 1] "Wetted perimeter"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Length dlengths[nFM] "Length of flow model"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Height[nFM + 1] roughnesses "Average height of surface asperities"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    input SI.Length[nFM] dheights
      "Height(states[2:nFM+1]) - Height(states[1:nFM])"
      annotation (Dialog(group="Inputs", tab="Internal Interface"));
    parameter Boolean allowFlowReversal = true
  "= true to allow flow reversal, false restricts to design direction (m_flows >= zeros(m))"
      annotation(Dialog(tab="Internal Interface", group="Advanced"), Evaluate=true);
    parameter SI.ReynoldsNumber Re_lam(max=Re_turb) = 2300
      "Laminar transition Reynolds number" annotation (Dialog(tab="Advanced"));
    parameter SI.ReynoldsNumber Re_turb(min=Re_lam) = 4000
      "Turbulent transition Reynolds number" annotation (Dialog(tab="Advanced"));
    parameter Boolean from_dp=momentumDynamics >= Modelica.Fluid.Types.Dynamics.SteadyStateInitial
      "= true, use m_flow = f(dp), otherwise dp = f(m_flow)"
      annotation (Dialog(tab="Advanced"), Evaluate=true);
    parameter SI.MassFlowRate m_flow_small=0.001 "Within regularization if |m_flows| < m_flow_small (may be wider for large discontinuities in static head)"
      annotation (Dialog(tab="Advanced",enable=not from_dp));
    parameter SI.AbsolutePressure dp_small = 1
      "Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)"
      annotation (Dialog(tab="Advanced",enable=from_dp));
  //   parameter Boolean use_Ts_film = false "=true for Ts_film = 0.5*(Ts_wall + Ts_fluid) else Ts_fluid" annotation(Dialog(tab="Advanced"));
    // Variables defined by model
    SI.MassFlowRate m_flows[nFM](start=m_flows_start, each stateSelect=if
          momentumDynamics == Dynamics.SteadyState then StateSelect.default else
          StateSelect.prefer,each min=if allowFlowReversal then -Modelica.Constants.inf else 0) "Mass flow rate across interfaces";
    // Total quantities
    SI.Momentum[nFM] Is "Momenta of flow segments";
    // Source/Sink terms
    SI.Force[nFM] Ibs "Flow of momentum across boundaries and source/sink in volumes";
    // Base Properties
    SI.Temperature[nFM + 1] Ts_fluid=Medium.temperature(states)
      "Fluid temperature";
  //   SI.Temperature[nFM + 1] Ts_film=if use_Ts_film then 0.5*(Ts_wall + Ts_fluid) else Ts_fluid "Film temperature";
  //   Medium.ThermodynamicState[nFM + 1] states_film=Medium.setState_pTX(
  //       Medium.pressure(states),
  //       Ts_film,
  //       Medium.X_default) "Film state";
    //Medium.X_default should be at leaste Medium.massFraction(state) but this doesn't seem to exist
    SI.ReynoldsNumber[nFM] Res "Reynolds number";
  //   SI.ReynoldsNumber[nFM] Res_film
  //     "Reynolds number with properties evaluated at film temperature";
  protected
    final parameter SI.ReynoldsNumber Re_center=0.5*(Re_lam + Re_turb)
      "Re smoothing transition center";
    final parameter SI.ReynoldsNumber Re_width=Re_turb - Re_center
      "Re smoothing transition width";
  initial equation
    if momentumDynamics == Dynamics.FixedInitial then
      m_flows = m_flows_start;
    elseif momentumDynamics == Dynamics.SteadyStateInitial then
      der(m_flows) = zeros(nFM);
    end if;
  equation
    for i in 1:nFM loop
      assert(m_flows[i] > -m_flow_small or allowFlowReversal, "Reverting flow occurs in flowModel even though allowFlowReversal is false");
    end for;
    // Total quantities
    Is = {m_flows[i]*dlengths[i] for i in 1:nFM};
    // Momentum balances
    if momentumDynamics == Dynamics.SteadyState then
      for i in 1:nFM loop
        0 = Ibs[i];
      end for;
    else
      for i in 1:nFM loop
        der(Is[i]) = Ibs[i];
      end for;
    end if;
    annotation (Documentation(info="<html>
</html>"),   Icon(graphics={Bitmap(extent={{-114,-100},{114,100}}, fileName="modelica://TRANSFORM/Resources/Images/Icons/FlowModel_dps.jpg")}));
  end PartialDistributedStaggeredFlow;

  partial model PartialMomentumBalance
    import Modelica.Fluid.Types.Dynamics;
    extends PartialDistributedStaggeredFlow;
    input Units.NonDim Ks_ab[nFM]=fill(0, nFM)
      "Minor loss coefficients. Flow in direction a -> b"
      annotation (Dialog(group="Inputs"));
    input Units.NonDim Ks_ba[nFM]=fill(0, nFM)
      "Minor loss coefficients. Flow in direction b -> a"
      annotation (Dialog(group="Inputs"));
  //   input SI.PressureDifference dps_add_ab[nFM]=fill(0, nFM)
  //     "Additional pressure losses. Flow in direction a -> b"
  //     annotation (Dialog(group="Inputs"));
  //   input SI.PressureDifference dps_add_ba[nFM]=fill(0, nFM)
  //     "Additional pressure losses. Flow in direction b -> a"
  //     annotation (Dialog(group="Inputs"));
    parameter Boolean use_I_flows=momentumDynamics <> Dynamics.SteadyState
      "= true to consider differences in flow of momentum through boundaries"
      annotation (Dialog(tab="Advanced"), Evaluate=true);
    parameter SI.Time taus[2]={0.01,0.01} "Time Constant for first order delay of {dps_K,dps_add}"
      annotation (Dialog(tab="Advanced"));
    // Source terms and forces to be defined by an extending model (zero if not used)
    SI.Force[nFM] I_flows "Flow of momentum across boundaries";
    SI.Force[nFM] Fs_p "Pressure forces";
    SI.Force[nFM] Fs_fg "Friction and gravity forces";
    SI.Pressure[nFM] dps_fg(start=dps_start) "Pressure drop between states";
    SI.PressureDifference dps_K[nFM] "Minor form-losses (K-loss)";
  //   SI.PressureDifference dps_add[nFM] "Minor additional pressure losses";
    Modelica.Blocks.Continuous.FirstOrder firstOrder_dps_K[nFM](
      each initType=Modelica.Blocks.Types.Init.InitialOutput,
      each y_start=0,
      each T=taus[1]) annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
  //   Modelica.Blocks.Continuous.FirstOrder firstOrder_dps_add[nFM](
  //     each initType=Modelica.Blocks.Types.Init.InitialOutput,
  //     each y_start=0,
  //     each T=taus[2])
  //     annotation (Placement(transformation(extent={{-80,-90},{-60,-70}})));
  equation
    Ibs = I_flows - Fs_p - Fs_fg;
    if use_I_flows then
      I_flows = {crossAreas[i]*Medium.density(states[i])*vs[i]*vs[i] - crossAreas[
        i + 1]*Medium.density(states[i + 1])*vs[i + 1]*vs[i + 1] for i in 1:nFM};
    else
      I_flows = zeros(nFM);
    end if;
    Fs_p = {0.5*(crossAreas[i] + crossAreas[i + 1])*(Medium.pressure(states[i+1]) -
      Medium.pressure(states[i])) for i in 1:nFM};
    Fs_fg = {(dps_fg[i] + firstOrder_dps_K[i].y)*0.5*(crossAreas[i] + crossAreas[i + 1]) for i in 1:nFM};
     dps_K = {noEvent(if m_flows[i] > 0 then 0.5*Ks_ab[i]*Medium.density(states[i])*vs[i]*vs[i] else
            - 0.5*Ks_ba[i]*Medium.density(states[i + 1])*vs[i + 1]*vs[i + 1]) for i in 1:nFM};
  //    dps_add = {noEvent(if m_flows[i] > 0 then dps_add_ab[i] else dps_add_ba[i]) for i in 1:nFM};
     firstOrder_dps_K[:].u = dps_K;
  //    firstOrder_dps_add[:].u = dps_add;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end PartialMomentumBalance;

  partial model PartialSinglePhase
    extends PartialMomentumBalance;
    TRANSFORM.Media.BaseProperties1Phase mediaProps[nFM + 1](redeclare package
        Medium = Medium, state=states) "Bulk fluid properties"
      annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
  //   TRANSFORM.Media.BaseProperties1Phase[nFM + 1] mediums_film(redeclare package
  //               Medium =
  //         Medium, state=states_film) "Film fluid properties"
  //     annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));
  equation
    for i in 1:nFM loop
      Res[i] =0.5 .* (dimensions[i] + dimensions[i + 1]) .* abs(m_flows[i]) ./ (
        0.25*(crossAreas[i] + crossAreas[i + 1]) .* (mediaProps[i].mu +
        mediaProps[i + 1].mu));
      //Res_film[i] = 0.5.*(dimensions[i] + dimensions[i+1]).*abs(m_flows[i])./(0.25*(crossAreas[i]+crossAreas[i+1]).*(mediums_film[i].mu+mediums_film[i+1].mu));
    end for;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end PartialSinglePhase;

  partial model PartialTwoPhase
    extends PartialMomentumBalance(
        redeclare replaceable package Medium =
          Modelica.Media.Water.StandardWater
        constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    replaceable model VoidFraction =
        TRANSFORM.Fluid.ClosureRelations.VoidFraction.Homogeneous_wSlipVelocity
                                                                                constrainedby
      TRANSFORM.Fluid.ClosureRelations.VoidFraction.PartialVoidFraction
      annotation (choicesAllMatching=true);
    TRANSFORM.Media.BaseProperties2Phase mediaProps[nFM + 1](redeclare package
        Medium = Medium, state=states,
      redeclare model VoidFraction = VoidFraction)
                                       "Bulk fluid properties"
      annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
  //   TRANSFORM.Media.BaseProperties2Phase[nFM+1] mediums_film(redeclare package
  //       Medium =
  //         Medium, state=states_film) "Film fluid properties"
  //     annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));
  equation
    for i in 1:nFM loop
      Res[i] =0.5 .* (dimensions[i] + dimensions[i + 1]) .* abs(m_flows[i]) ./ (
        0.25*(crossAreas[i] + crossAreas[i + 1]) .* (mediaProps[i].mu +
        mediaProps[i + 1].mu));
      //Res_film[i] = 0.5.*(dimensions[i] + dimensions[i+1]).*abs(m_flows[i])./(0.25*(crossAreas[i]+crossAreas[i+1]).*(mediums_film[i].mu+mediums_film[i+1].mu));
    end for;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end PartialTwoPhase;

  model SinglePhase_Developed_2Region_NumStable
    "Single Phase | Fully Developed | 2 Region - Laminar & Turbulent | Numerically Stable Method"
    extends PartialSinglePhase;
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_IN_con[
      nFM] IN_con(
      length=dlengths,
      diameter_a=dimensions[1:nFM],
      diameter_b=dimensions[2:nFM + 1],
      crossArea_a=crossAreas[1:nFM],
      crossArea_b=crossAreas[2:nFM + 1],
      roughness_a=roughnesses[1:nFM],
      roughness_b=roughnesses[2:nFM + 1],
      each Re_turbulent=Re_turb)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_IN_var[
      nFM] IN_var(
      rho_a=ds_a,
      rho_b=ds_b,
      mu_a=mus_a,
      mu_b=mus_b)
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
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
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_IN_var
      IN_var_nominal(
      rho_a=d_nominal,
      rho_b=d_nominal,
      mu_a=mu_nominal,
      mu_b=mu_nominal);
    SI.AbsolutePressure dp_fric_nominal=sum(
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_DP(
        IN_con,
        IN_var_nominal,
        m_flow_nominal,
        m_flow_small)) "pressure loss for nominal conditions";
  initial equation
    // dp_nominal = dp_fric_nominal + g_n*sum(dheights)*d_nominal;
      dp_nominal = dp_fric_nominal + sum(g_eff * dheights)*d_nominal;
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_MFLOW(
          IN_con,
          IN_var,
          dps_fg - {g_eff[i]*dheights[i]*ds_a[i] for i in 1:nFM},
          dp_small/(nFM)), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_eff.*dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_DP(
          IN_con,
          IN_var,
          m_flows,
          m_flow_small) + {g_eff[i]*dheights[i]*ds_a[i] for i in 1:nFM}, simplified=dp_nominal/m_flow_nominal*m_flows + g_eff.*dheights*d_nominal);
      end if;
    else
      // regularization for discontinuous flow reversal and static head
      if from_dp then
        // and not dp_is_zero then
        m_flows =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_MFLOW_staticHead(
          IN_con,
          IN_var,
          dps_fg,
          dp_small/(nFM),g_eff.*dheights), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_eff.*dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_DP_staticHead(
          IN_con,
          IN_var,
          m_flows,
          m_flow_small,
          g_eff.*dheights), simplified=dp_nominal/m_flow_nominal*m_flows + g_eff.*dheights*d_nominal);
      end if;
    end if;
      annotation (Dialog(tab="Advanced", group="Nominal Conditions",enable=false),
                Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end SinglePhase_Developed_2Region_NumStable;

  model SinglePhase_Turbulent_MSL
    "Single Phase | Turbulent | Numerically Stable Method"
    extends PartialSinglePhase;
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_IN_con[
      nFM] IN_con(
      length=dlengths,
      diameter_a=dimensions[1:nFM],
      diameter_b=dimensions[2:nFM + 1],
      crossArea_a=crossAreas[1:nFM],
      crossArea_b=crossAreas[2:nFM + 1],
      roughness_a=roughnesses[1:nFM],
      roughness_b=roughnesses[2:nFM + 1],
      each Re_turbulent=Re_turb)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_IN_var[
      nFM] IN_var(
      rho_a=ds_a,
      rho_b=ds_b,
      mu_a=mus_a,
      mu_b=mus_b)
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
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
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_IN_var
      IN_var_nominal(
      rho_a=d_nominal,
      rho_b=d_nominal,
      mu_a=mu_nominal,
      mu_b=mu_nominal);
    SI.AbsolutePressure dp_fric_nominal=sum(
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_DP(
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_MFLOW(
          IN_con,
          IN_var,
          dps_fg - {g_n*dheights[i]*ds_a[i] for i in 1:nFM},
          dp_small/(nFM)), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
          dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_DP(
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_MFLOW_staticHead(
          IN_con,
          IN_var,
          dps_fg,
          dp_small/(nFM),
          g_n*dheights), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
          dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarAndQuadraticTurbulent_MSL.dp_DP_staticHead(
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
  end SinglePhase_Turbulent_MSL;

  model SinglePhase_Developed_2Region_Simple
    "Single Phase | Fully Developed | 2 Region - Laminar & Turbulent | Simple Method"
    import TRANSFORM;
    extends PartialSinglePhase;
    input TRANSFORM.Units.NonDim[nFM] fRe2=
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.fRe2_SinglePhase_2Region(
        Res,
        dimensionsAvg,
        roughnessesAvg,
        Re_center,
        Re_width) "Turbulent Friction factor"
      annotation (Dialog(group="Inputs"));
    SI.Pressure[nFM] dps_f "Frictional pressure drop";
    SI.Pressure[nFM] dps_g "Gravitational pressure drop";
  protected
    SI.DynamicViscosity mus[nFM];
    SI.DynamicViscosity mus_a[nFM]=mediaProps[1:nFM].mu;
    SI.DynamicViscosity mus_b[nFM]=mediaProps[2:nFM + 1].mu;
    SI.Density ds[nFM];
    SI.Density ds_a[nFM]=mediaProps[1:nFM].d;
    SI.Density ds_b[nFM]=mediaProps[2:nFM + 1].d;
    SI.Length dimensionsAvg[nFM];
    SI.Length roughnessesAvg[nFM];
  equation
    dps_fg = dps_f+dps_g;
    for i in 1:nFM loop
      dps_f[i] = 0.5*fRe2[i]*dlengths[i]*mus[i]^2/(dimensionsAvg[i]*dimensionsAvg[
        i]*dimensionsAvg[i]*ds[i])*noEvent(if m_flows[i] >= 0 then +1 else -1);
      dps_g[i] = g_n*dheights[i]*0.5*(Medium.density(states[i]) + Medium.density(
        states[i + 1]));
    end for;
    for i in 1:nFM loop
      dimensionsAvg[i] = 0.5*(dimensions[i] + dimensions[i + 1]);
      roughnessesAvg[i] = 0.5*(roughnesses[i] + roughnesses[i + 1]);
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
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end SinglePhase_Developed_2Region_Simple;

  model Nominal_Linear "Set nominal dp and m_flow"
    extends PartialSinglePhase;
    parameter SI.AbsolutePressure dp_nominal=1 "Nominal pressure loss across pipe";
    parameter SI.MassFlowRate m_flow_nominal=1 "Mass flow rate for dp_nominal";
    SI.Pressure[nFM] dps_f "Frictional pressure drop";
    SI.Pressure[nFM] dps_g "Gravitational pressure drop";
  equation
    dps_fg = dps_f + dps_g;
    for i in 1:nFM loop
      dps_f[i] =dp_nominal/m_flow_nominal/nFM*m_flows[i];
      dps_g[i] =g_n*dheights[i]*0.5*(Medium.density(states[i]) + Medium.density(
        states[i + 1]));
    end for;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Nominal_Linear;

  model TwoPhase_Developed_2Region_NumStable
    "Two Phase | Fully Developed | 2 Region - Laminar & Turbulent | Numerically Stable Method"
    import TRANSFORM;
    extends PartialTwoPhase;
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_IN_con[
      nFM] IN_con(
      length=dlengths,
      diameter_a=dimensions[1:nFM],
      diameter_b=dimensions[2:nFM + 1],
      crossArea_a=crossAreas[1:nFM],
      crossArea_b=crossAreas[2:nFM + 1],
      roughness_a=roughnesses[1:nFM],
      roughness_b=roughnesses[2:nFM + 1],
      each Re_turbulent=Re_turb)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_IN_var[
      nFM] IN_var(
      rho_a=ds_a,
      rho_b=ds_b,
      mu_a=mus_a,
      mu_b=mus_b)
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
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
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_IN_var
      IN_var_nominal(
      rho_a=d_nominal,
      rho_b=d_nominal,
      mu_a=mu_nominal,
      mu_b=mu_nominal);
    SI.AbsolutePressure dp_fric_nominal=sum(
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_DP(
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_MFLOW(
          IN_con,
          IN_var,
          dps_fg - {g_n*dheights[i]*ds_a[i] for i in 1:nFM},
          dp_small/(nFM)), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
          dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_DP(
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_MFLOW_staticHead(
          IN_con,
          IN_var,
          dps_fg,
          dp_small/(nFM),
          g_n*dheights), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
          dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.SinglePhase.LaminarTurbulent_MSLDetailed.dp_DP_staticHead(
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
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>Same as SinglePhase_Developed_2Region_NumStable but extends from PartialTwoPhase which uses alphaV for calculation of density.</p>
</html>"));
  end TwoPhase_Developed_2Region_NumStable;

  model TwoPhase_Developed_2Region_NumStable_alternate
    "Two Phase | Fully Developed | 2 Region - Laminar & Turbulent | Numerically Stable Method"
    extends PartialTwoPhase;
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_IN_con[
      nFM] IN_con(
      length=dlengths,
      diameter_a=dimensions[1:nFM],
      diameter_b=dimensions[2:nFM + 1],
      crossArea_a=crossAreas[1:nFM],
      crossArea_b=crossAreas[2:nFM + 1],
      roughness_a=roughnesses[1:nFM],
      roughness_b=roughnesses[2:nFM + 1],
      each Re_turbulent=Re_turb)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_IN_var[
      nFM] IN_var(
      rho_a=ds_a,
      rho_b=ds_b,
      rho_lsat_a=ds_lsat_a,
      rho_lsat_b=ds_lsat_b,
      rho_vsat_a=ds_vsat_a,
      rho_vsat_b=ds_vsat_b,
      mu_a=mus_a,
      mu_b=mus_b,
      mu_lsat_a=mus_lsat_a,
      mu_lsat_b=mus_lsat_b,
      mu_vsat_a=mus_vsat_a,
      mu_vsat_b=mus_vsat_b,
      x_abs_a=x_abs_a,
      x_abs_b=x_abs_b)
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
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
    SI.Density ds_lsat_a[nFM]=if use_d_nominal then fill(d_nominal, nFM) else
        mediaProps[1:nFM].rho_lsat;
    SI.Density ds_lsat_b[nFM]=if use_d_nominal then fill(d_nominal, nFM) else
        mediaProps[2:nFM + 1].rho_lsat;
    SI.Density ds_vsat_a[nFM]=if use_d_nominal then fill(d_nominal, nFM) else
        mediaProps[1:nFM].rho_vsat;
    SI.Density ds_vsat_b[nFM]=if use_d_nominal then fill(d_nominal, nFM) else
        mediaProps[2:nFM + 1].rho_vsat;
    SI.DynamicViscosity mus_lsat_a[nFM]=if use_mu_nominal then fill(mu_nominal, nFM)
         else mediaProps[1:nFM].mu_lsat;
    SI.DynamicViscosity mus_lsat_b[nFM]=if use_mu_nominal then fill(mu_nominal, nFM)
         else mediaProps[2:nFM + 1].mu_lsat;
    SI.DynamicViscosity mus_vsat_a[nFM]=if use_mu_nominal then fill(mu_nominal, nFM)
         else mediaProps[1:nFM].mu_vsat;
    SI.DynamicViscosity mus_vsat_b[nFM]=if use_mu_nominal then fill(mu_nominal, nFM)
         else mediaProps[2:nFM + 1].mu_vsat;
    SIadd.NonDim x_abs_a[nFM]= mediaProps[1:nFM].x_abs;
    SIadd.NonDim x_abs_b[nFM]= mediaProps[2:nFM + 1].x_abs;
  protected
    TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_IN_var
      IN_var_nominal(
      rho_a=d_nominal,
      rho_b=d_nominal,
      rho_lsat_a=d_nominal,
      rho_lsat_b=d_nominal,
      rho_vsat_a=d_nominal,
      rho_vsat_b=d_nominal,
      mu_a=mu_nominal,
      mu_b=mu_nominal,
      mu_lsat_a=mu_nominal,
      mu_lsat_b=mu_nominal,
      mu_vsat_a=mu_nominal,
      mu_vsat_b=mu_nominal,
      x_abs_a=0,
      x_abs_b=0);
    SI.AbsolutePressure dp_fric_nominal=sum(
        TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_DP(
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_MFLOW(
          IN_con,
          IN_var,
          dps_fg - {g_n*dheights[i]*ds_a[i] for i in 1:nFM},
          dp_small/(nFM)), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
          dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_DP(
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
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_MFLOW_staticHead(
          IN_con,
          IN_var,
          dps_fg,
          dp_small/(nFM),
          g_n*dheights), simplified=m_flow_nominal/dp_nominal*(dps_fg - g_n*
          dheights*d_nominal));
      else
        dps_fg =homotopy(actual=
          TRANSFORM.Fluid.ClosureRelations.PressureLoss.Functions.TubesAndConduits.TwoPhase.LaminarTurbulent_MSLDetailed.dp_DP_staticHead(
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
  end TwoPhase_Developed_2Region_NumStable_alternate;

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
end DistributedPipe_1D_ver1;
