within TRANSFORM.Math.Easing.Elastic;
function easeIn "Elastic | Ease In"
  extends PartialEasing;

protected
  Real scaledX =  x/deltax;
  Real y_int;

  Real a = 2*Modelica.Constants.pi/3;

algorithm
  if scaledX <= -0.999999999 then
    y_int := 0;
  elseif scaledX >= 0.999999999 then
    y_int := 1;
  else
    y_int := -2^(10 * scaledX - 10) * sin((scaledX * 10 - 10.75) * a);
  end if;
  y := pos*y_int + (1 - y_int)*neg;

  annotation (smoothOrder=1);
end easeIn;
