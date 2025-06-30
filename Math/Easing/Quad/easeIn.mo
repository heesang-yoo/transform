within TRANSFORM.Math.Easing.Quad;
function easeIn "Quad | Ease In"
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
    y_int := (scaledX+1)^2/4;
  end if;
  y := pos*y_int + (1 - y_int)*neg;

  annotation (smoothOrder=1);
end easeIn;
