within TRANSFORM.Fluid.Pipes.Examples;
model pipeflow_1 "Comparing a circular with a non-circular pipe"
  import TRANSFORM;
  extends TRANSFORM.Icons.Example;
  replaceable package Medium = Modelica.Media.Water.StandardWater;
  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
                                     annotation(Placement(transformation(origin = {-80, 80}, extent = {{-10, -10}, {10, 10}})));
  TRANSFORM.Fluid.BoundaryConditions.Boundary_pT boundary(
    redeclare package Medium = Medium,
    nPorts=1,
    p=10.0e5,
    T=293.15) annotation (Placement(transformation(origin={-34,36},extent={{-10,
            -10},{10,10}},
        rotation=270)));
  TRANSFORM.Fluid.BoundaryConditions.MassFlowSource_T massflowsink1(
    redeclare package Medium = Medium,
    m_flow=0.0,
    T=293.15,
    nPorts=1) annotation (Placement(transformation(origin={-34,-54},
                                                                   extent={{10,
            -10},{-10,10}},
        rotation=270)));
TRANSFORM.Fluid.Pipes.GenericPipe_MultiTransferSurface
            circular_pipe(
  redeclare package Medium = Medium,
  energyDynamics=system.energyDynamics,
  m_flow_a_start=system.m_flow_start,
  momentumDynamics=system.momentumDynamics,
  exposeState_b=true,
  exposeState_a=false,
    p_a_start=1000000,
    p_b_start=1000000,
    T_a_start=293.15,
    redeclare model Geometry =
        TRANSFORM.Fluid.ClosureRelations.Geometry.Models.DistributedVolume_1D.StraightPipe
        (
        dimension=0.01,
        length=100,
        angle=1.5707963267949))      annotation (Placement(transformation(
        origin={-34,-12},extent={{-10,-10},{10,10}},
        rotation=270)));
equation
  connect(boundary.ports[1], circular_pipe.port_a) annotation(Line(points={{-34,26},
          {-34,-2}},                                                                                                    color = {0, 127, 255}));
  connect(circular_pipe.port_b, massflowsink1.ports[1])
    annotation (Line(points={{-34,-22},{-34,-44}}, color={0,127,255}));
  annotation(experiment(StopTime=1),
  Documentation(info="<html>
  <h4> simple pipe flow with gravity </h4>
<p>
In this example two pipes are used to demonstrate the use of circular (default) and non-circular pipes,
where the topmost pipe is circular with a length of 100 m and an inner diameter of 10 mm and the second pipe
is a circular ring pipe with inner diameter of 5 mm and an outer diameter of 15 mm.
</p>
<p>
Both pipes are connected to a pT source (water, 293.15 K, 10 bar) and a mass flow sink (0.1 kg/s inflow).
</p>
<p>
Although the hydraulic diameter of both pipes are the same, the different cross sections lead to different
velocities and by this different outlet pressures (7.324 bar for the circular pipe versus 9.231 bar for the
circular ring pipe).
</p>
</html>", revisions="<html>
<ul>
<li>
January 6, 2015 by Alexander T&auml;schner:<br/>
First implementation.
</li>
</ul>
</html>"));
end pipeflow_1;
