within TRANSFORM.Fluid.Pipes;
model GenericPipe_MultiTransferSurface_test1
  extends GenericPipe_MultiTransferSurface(
    redeclare input SI.Acceleration g_n = g_ext[3]);

  input SI.Acceleration g_ext[3] = {0,0,-Modelica.Constants.g_n};
  annotation (
    defaultComponentName="pipe",
    Documentation(info="<html>
    <h4>extended from GenericPipe_MultiTransferSurface</h4>
    </html>"));
end GenericPipe_MultiTransferSurface_test1;
