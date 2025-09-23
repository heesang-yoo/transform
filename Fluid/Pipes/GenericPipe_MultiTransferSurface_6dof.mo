within TRANSFORM.Fluid.Pipes;
model GenericPipe_MultiTransferSurface_6dof
  extends GenericPipe_MultiTransferSurface(
    redeclare input SI.Acceleration g_n = g_ext[3]); // 바인딩 override

  input SI.Acceleration g_ext[3] = {0,0,-Modelica.Constants.g_n};

end GenericPipe_MultiTransferSurface_6dof;
