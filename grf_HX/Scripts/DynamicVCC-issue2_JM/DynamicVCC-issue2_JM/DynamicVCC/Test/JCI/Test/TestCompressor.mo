within DynamicVCC.Test.JCI.Test;
model TestCompressor
  extends Modelica.Icons.Example;

  inner DynamicVCC.Components.System system(
  redeclare package Medium=Medium);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  DynamicVCC.Test.JCI.Compressor compressor(
  redeclare package Medium=Medium,
  Vs=2.9438902e-05,
  m_flow_nominal=0.0279,
  m_flow_init=0,
  speed_nominal=42);

  Modelica.Fluid.Sources.Boundary_ph suction(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true,
  use_h_in=true);

  Modelica.Fluid.Sources.Boundary_ph discharge(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true);

equation

  connect(suction.ports[1],compressor.port_a);
  connect(compressor.port_b,discharge.ports[1]);

  suction.p_in=10e5;
  discharge.p_in=20e5;
  suction.h_in=4.3e5;
  compressor.speed=58;
  compressor.T_amb=310.5;

  annotation (experiment(
      StartTime=1,
      StopTime=100,
      __Dymola_Algorithm="Dassl"));
end TestCompressor;
