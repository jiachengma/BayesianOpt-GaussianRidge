within DynamicVCC.Test.JCI.Test.TPWL;
model LinearizeCompressor
  extends Modelica.Icons.Example;

  input Real u[4](start={0.0494,4.3762e5,23.5769e5,58.3});
  output Real y[4];

  inner DynamicVCC.Components.System system(
  redeclare package Medium=Medium);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  DynamicVCC.Test.JCI.Compressor compressor(
  redeclare package Medium=Medium,
  Vs=2.9438902e-05,
  m_flow_nominal=0.051,
  m_flow_init=0.05,
  speed_nominal=58);

  Modelica.Fluid.Sources.MassFlowSource_h suction(
  redeclare package Medium=Medium,
  nPorts=1,
  use_m_flow_in=true,
  use_h_in=true);

  Modelica.Fluid.Sources.Boundary_ph discharge(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true);

equation

  connect(suction.ports[1],compressor.port_a);
  connect(compressor.port_b,discharge.ports[1]);

  suction.m_flow_in=u[1];
  suction.h_in=u[2];
  discharge.p_in=u[3];
  compressor.speed=u[4];
  compressor.T_amb=305;

  y[1]=compressor.p_suc;
  y[2]=compressor.h_dis;
  y[3]=compressor.Pwr;
  y[4]=compressor.m_flow;

  annotation (experiment(
      StartTime=1,
      StopTime=100,
      __Dymola_Algorithm="Dassl"));
end LinearizeCompressor;
