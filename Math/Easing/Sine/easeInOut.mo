within TRANSFORM.Math.Easing.Sine;
function easeInOut "Sine | Ease In & Out"
  extends PartialEasing;

protected
  Real scaledX =  x/deltax;
  Real y_int;

algorithm
  if scaledX <= -0.999999999 then
    y_int := 0;
  elseif scaledX >= 0.999999999 then
    y_int := 1;
  else
    y_int := -cos(0.5*(scaledX+1)*Modelica.Constants.pi)/2 + 0.5;
  end if;
  y := pos*y_int + (1 - y_int)*neg;

  annotation (smoothOrder=1);
end easeInOut;
