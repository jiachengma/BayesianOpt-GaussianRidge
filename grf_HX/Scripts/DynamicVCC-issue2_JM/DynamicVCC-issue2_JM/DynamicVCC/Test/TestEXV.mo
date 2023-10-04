within DynamicVCC.Test;
model TestEXV
  extends Modelica.Icons.Example;

  inner DynamicVCC.Components.System system(
  m_flow_init=m_flow_init,
  m_flow_nominal=m_flow_nominal,
  enableReverseFlow=true);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  parameter Medium.MassFlowRate m_flow_init=0.04765;
  parameter Medium.MassFlowRate m_flow_nominal=m_flow_init;

  DynamicVCC.Examples.GreenspeedASHP.EXV exv(
    redeclare package Medium = Medium,
    Av=3.14e-6,
    final dp_nominal=10e5,
    opening_init=0.16);

  Modelica.Fluid.Sources.Boundary_ph source(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true,
  use_h_in=true);

  Modelica.Fluid.Sources.Boundary_ph sink(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true);

equation
  connect(source.ports[1],exv.port_a);
  connect(exv.port_b,sink.ports[1]);

  source.p_in=3392977;
  source.h_in=243818;
  sink.p_in=1226031;

  exv.opening=0.15;

  annotation (experiment(
      StartTime=1,
      StopTime=100,
      Interval=10,
      __Dymola_Algorithm="Dassl"));
end TestEXV;
