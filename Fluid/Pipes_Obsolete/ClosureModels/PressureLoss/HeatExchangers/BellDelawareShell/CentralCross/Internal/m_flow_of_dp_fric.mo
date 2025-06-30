within TRANSFORM.Fluid.Pipes_Obsolete.ClosureModels.PressureLoss.HeatExchangers.BellDelawareShell.CentralCross.Internal;
function m_flow_of_dp_fric
  "Calculate mass flow rate as a function of pressure drop"
  extends Modelica.Icons.Function;
  //input records
  input
    TRANSFORM.Fluid.Pipes_Obsolete.ClosureModels.PressureLoss.HeatExchangers.BellDelawareShell.CentralCross.dp_IN_con
    IN_con "Input record for function dp_overall_MFLOW"
    annotation (Dialog(group="Constant inputs"));
  input
    TRANSFORM.Fluid.Pipes_Obsolete.ClosureModels.PressureLoss.HeatExchangers.BellDelawareShell.CentralCross.dp_IN_var
    IN_var "Input record for function dp_overall_MFLOW"
    annotation (Dialog(group="Variable inputs"));
  input SI.Pressure dp_fric
    "Pressure loss due to friction (dp = port_a.p - port_b.p)";
  //Outputs
  output SI.MassFlowRate m_flow;
  output Real dm_flow_ddp_fric "Derivative of mass flow rate with dp_fric";
protected
  Real a = IN_con.s1/IN_con.d_B;
  Real b = IN_con.s2/IN_con.d_B;
  Real c = ((a/2)^2 + b^2)^(0.5);
  SI.Length e = (if IN_con.toggleStaggered then
                  (if b >= 0.5*(2*a+1)^(0.5) then (a - 1)*IN_con.d_B else (c - 1)*IN_con.d_B)
                 else (a - 1)*IN_con.d_B);
  SI.Length L_E = 2*IN_con.e1 + e*IN_con.nes;
  SI.Area A_E = IN_con.S*L_E;
  Real gamma = 2*Modelica.Math.acos(1 - 2*IN_con.H/IN_con.D_l)*180/pi;
  Real f_L=
    Internal.f_L_baffact(
      gamma,
      IN_con.H,
      IN_con.D_l,
      IN_con.D_i,
      IN_con.d_B,
      IN_con.d_B,
      IN_con.n_T,
      IN_con.n_W,
      A_E);
  SI.Area A_B = if e<(IN_con.D_i-IN_con.DB) then IN_con.S*(IN_con.D_i-IN_con.DB-e) else 0;
  Real R_B = A_B/A_E;
  Real R_S = IN_con.n_s/IN_con.n_MR;
  Real beta;
  Real f_B;
  Real f_zL;
  Real f_zt;
  Real df_zL_dm_flow;
  Real df_zt_dm_flow;
  Real epsilon;
  Real depsilon_dm_flow;
  Real error = 1;
  Real errTol = 1e-3;
  Real iter = 0;
  Real itMax = 10;
  SI.ReynoldsNumber Re_old = 1e4;
  SI.DynamicViscosity mu "Upstream viscosity";
  SI.Density rho "Upstream density";
  Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
  SI.ReynoldsNumber Re "Reynolds number";
  Real dRe_ddp "dRe/ddp";
  Real aux1;
  Real aux2;
  Real aux3;
algorithm
  // Determine upstream density and upstream viscosity
  rho := if dp_fric >= 0 then IN_var.rho_a else IN_var.rho_b;
  mu  := if dp_fric >= 0 then IN_var.mu_a else IN_var.mu_b;
  aux3 :=IN_con.d_B/(mu*A_E);
  aux1:= 2*rho*IN_con.d_B*IN_con.d_B/(mu*mu);
  lambda2 := abs(dp_fric)*aux1*IN_con.nNodes;
  while error > errTol and iter < itMax loop
    iter := iter + 1;
    beta :=if Re_old < 100 then 4.5 else 3.7;
    if R_S == 0 then
      f_B :=exp(-beta*R_B);
    elseif R_S < 0.5 then
      f_B := exp(-beta*R_B*(1 - (2*R_S)^(1/3)));
    else
      f_B := 1;
    end if;
    (f_zL,f_zt,df_zL_dm_flow,df_zt_dm_flow) :=
      Internal.f_LeakFactors(
        Re_old,
        a,
        b,
        mu,
        IN_var.mu_w);
    (epsilon) :=
      Internal.DragCoeff(
        IN_con.toggleStaggered,
        Re_old,
        a,
        b,
        c,
        f_zL,
        f_zt,
        df_zL_dm_flow,
        df_zt_dm_flow,
        aux3);
    aux2 :=epsilon*IN_con.n_MR*f_L*f_B;
    // Determine Re
    Re := sqrt(lambda2/aux2);
    error := abs(Re - Re_old);
    Re_old := Re;
  end while;
  (epsilon) :=
    Internal.DragCoeff(
      IN_con.toggleStaggered,
      Re,
      a,
      b,
      c,
      f_zL,
      f_zt,
      df_zL_dm_flow,
      df_zt_dm_flow,
      aux3);
  aux2 :=epsilon*IN_con.n_MR*f_L*f_B;
  dRe_ddp := 0.5*aux1*IN_con.nNodes/aux2/sqrt(lambda2/aux2);
  // Determine mass flow rate
  m_flow := mu*A_E/IN_con.d_B*(if dp_fric >= 0 then Re else -Re);
  // Determine derivative of mass flow rate with dp_fric
  dm_flow_ddp_fric := mu*A_E/IN_con.d_B*dRe_ddp;
  annotation(smoothOrder=1);
end m_flow_of_dp_fric;
