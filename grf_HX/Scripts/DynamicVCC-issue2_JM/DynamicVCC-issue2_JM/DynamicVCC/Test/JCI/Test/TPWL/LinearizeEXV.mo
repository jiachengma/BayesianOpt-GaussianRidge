within DynamicVCC.Test.JCI.Test.TPWL;
model LinearizeEXV
  extends Modelica.Icons.Example;

  input Real u[4](start={0.0494,2.53e5,9.84e5,0.45});

  output Real y[3];

  inner DynamicVCC.Components.System system(
  m_flow_init=0.05,
  m_flow_nominal=0.05,
  massDynamics=DynamicVCC.Components.Types.Dynamics.DynamicFree_init,
  energyDynamics=DynamicVCC.Components.Types.Dynamics.DynamicFree_init,
  momentumDynamics=DynamicVCC.Components.Types.Dynamics.SteadyState,
  enableReverseFlow=true);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  DynamicVCC.Test.JCI.EXV exv(
  redeclare package Medium=Medium,
  Av=3.1416e-6,
  final dp_nominal=10e5);

  Modelica.Fluid.Sources.MassFlowSource_h source(
  redeclare package Medium=Medium,
  nPorts=1,
  use_m_flow_in=true,
  use_h_in=true);

  Modelica.Fluid.Sources.Boundary_ph sink(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true);

equation
  connect(source.ports[1],exv.port_a);
  connect(exv.port_b,sink.ports[1]);

  source.m_flow_in=u[1];
  source.h_in=u[2];
  sink.p_in=u[3];
  exv.opening=u[4];

  y[1]=exv.port_a.p;
  y[2]=exv.port_b.h_outflow;
  y[3]=exv.m_flow;

end LinearizeEXV;
