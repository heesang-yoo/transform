within TRANSFORM.Math;
function radialspace "Equal volume radial spacing"
  input Real r_inner "Inner radius";
  input Real r_outer "Outer radius";
  input Integer n "# of spacings";
  output Real rs[n] "points";
algorithm
  rs[1] :=r_inner;
  for i in 2:n loop
    rs[i] :=sqrt((r_outer^2 - r_inner^2)/(n-1) + rs[i - 1]^2);
  end for;
  annotation (Documentation(info="<html>
<p><span style=\"font-family: Courier New;\">For example:</span></p>
<p><span style=\"font-family: Courier New;\">radialspace(0,1,5) = {0,0.5,0.707,0.866,1}</span></p>
</html>"));
end radialspace;
