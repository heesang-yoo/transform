within TRANSFORM.Fluid;
package Sensors
  extends TRANSFORM.Icons.SensorsPackage;

  annotation (preferredView="info", Documentation(info="<html>
<p align = justify>
Package <b>Sensors</b> consists of idealized sensor components that
provide variables of a medium model and/or fluid ports as
output signals. These signals can be, e.g., further processed
with components of the Modelica.Blocks library.
Also more realistic sensor models can be built, by further
processing (e.g., by attaching block Modelica.Blocks.FirstOrder to
model the time constant of the sensor).
</p>

<p align = justify>For the thermodynamic state variables temperature, specific enthalpy, specific entropy and density
the fluid library provides two different types of sensors: <b>regular one port</b> and <b>two port</b> sensors.</p>

<ul>
<li>The <b>regular one port</b> sensors have the advantage of easy introduction and removal from a model, as no connections have to be broken.
A potential drawback is that the obtained value jumps as flow reverts.
</li>

<li>The <b>two port</b> sensors offer the advantages of an adjustable regularized step function around zero flow.
Moreover the obtained result is restricted to the value flowing into port_a if allowFlowReversal is false.</li>
</ul>

<p>
<a href=\"modelica://Modelica.Fluid.Examples.Explanatory.MeasuringTemperature\">Modelica.Fluid.Examples.Explanatory.MeasuringTemperature</a>
demonstrates the differences between one- and two-port sensor at hand of a
simple example.
</p>
</html>", revisions="<html>
<ul>
<li><i>22 Dec 2008</i>
    by R;uumldiger Franke<br>
    <ul>
    <li>flow sensors based on Interfaces.PartialTwoPort</li>
    <li>adapted docu to stream connectors, i.e., less need for two port sensors</li>
    </ul>
<li><i>4 Dec 2008</i>
    by Michael Wetter<br>
       included sensors for trace substance</li>
<li><i>31 Oct 2007</i>
    by Carsten Heinrich<br>
       updated sensor models, included one and two port sensors for thermodynamic state variables</li>
</ul>
</html>"));
end Sensors;
