within DynamicVCC.Test;
model TestCompressor
  extends Modelica.Icons.Example;

  inner DynamicVCC.Components.System system(
  redeclare package Medium=Medium,
  m_flow_nominal=m_flow_nominal);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  parameter Medium.MassFlowRate m_flow_init=0.048;
  parameter Medium.MassFlowRate m_flow_nominal=m_flow_init;
  parameter SI.Frequency speed_nominal=52;

  DynamicVCC.Examples.GreenspeedASHP.Compressor compressor(
    redeclare final package Medium = Medium,
    Vs=2.3765e-5,
    speed_nominal=speed_nominal,
    UA=16.92);

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

  suction.p_in=1167704;
  discharge.p_in=3415010;
  suction.h_in=424893;
  compressor.speed=51.75;
  compressor.T_amb=286.6;

  annotation (experiment(
      StartTime=1,
      StopTime=100,
      __Dymola_Algorithm="Dassl"));
end TestCompressor;
