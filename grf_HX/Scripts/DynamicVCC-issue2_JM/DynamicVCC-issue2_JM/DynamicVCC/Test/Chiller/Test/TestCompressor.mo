within DynamicVCC.Test.Chiller.Test;
model TestCompressor
  extends Modelica.Icons.Example;

  inner DynamicVCC.Components.System system(
  p_max=12e5,
  p_min=2e5,
  h_max=5e5,
  h_min=1e5,
  T_max=340,
  T_min=260,
  EnableReverseFlow=false);

  replaceable package Medium=Modelica.Media.R134a.R134a_ph;

  DynamicVCC.Test.Chiller.Compressor compressor(redeclare package Medium=Medium);

  // Suction
  Modelica.Fluid.Sources.Boundary_ph source(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true,
  use_h_in=true);

  // Discharge
  Modelica.Fluid.Sources.Boundary_ph sink(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true);

  Modelica.Blocks.Sources.CombiTimeTable BC_Cond(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Cond",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Cond.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable BC_Evap(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Evap",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Evap.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:9);

equation

  compressor.gamma=1;

  source.p_in=Mea.y[1];
  source.h_in=Mea.y[5];
  sink.p_in=Mea.y[2];

  connect(source.ports[1],compressor.port_a);
  connect(compressor.port_b,sink.ports[1]);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=2000,
      StopTime=5000,
      __Dymola_Algorithm="Dassl"));
end TestCompressor;
