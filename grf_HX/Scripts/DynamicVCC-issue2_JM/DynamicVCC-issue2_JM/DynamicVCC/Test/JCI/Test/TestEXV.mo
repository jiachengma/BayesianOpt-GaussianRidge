within DynamicVCC.Test.JCI.Test;
model TestEXV
  extends Modelica.Icons.Example;

   inner DynamicVCC.Components.System system(
  m_flow_init=0,
  m_flow_nominal=0.051,
  massDynamics=DynamicVCC.Components.Types.Dynamics.DynamicFree_init,
  energyDynamics=DynamicVCC.Components.Types.Dynamics.DynamicFree_init,
  momentumDynamics=DynamicVCC.Components.Types.Dynamics.SteadyState,
  enableReverseFlow=true);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  DynamicVCC.Test.JCI.EXV exv(
  redeclare package Medium=Medium,
  Av=3.1416e-6,
  final dp_nominal=10e5);

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

  source.p_in=20e5;
  source.h_in=2.5e5;
  sink.p_in=20e5;

  exv.opening=0.5;

end TestEXV;
