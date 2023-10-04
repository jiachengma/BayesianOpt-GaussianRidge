within DynamicVCC.Test.Chiller.Test;
model TestTXV
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

  DynamicVCC.Test.Chiller.Bulb bulb(C=50);

  DynamicVCC.Test.Chiller.TXV txv(
  redeclare package Medium=Medium,
  dp_nominal=5e5,
  m_flow_init=2,
  m_flow_nominal=2);

  Medium.ThermodynamicState state_suc; //suction state

  //Source
  Modelica.Fluid.Sources.Boundary_ph source(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true,
  use_h_in=true);

  //Sink
  Modelica.Fluid.Sources.Boundary_ph sink(
  redeclare package Medium=Medium,
  nPorts=1,
  use_p_in=true);

  Modelica.Blocks.Sources.CombiTimeTable BC_Cond(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Cond",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Cond.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable BC_Evap(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Evap",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Evap.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:9);

equation
  state_suc=Medium.setState_ph(Mea.y[1],Mea.y[5]);
  bulb.u=Medium.temperature(state_suc);
  txv.Q_dot_m=0;

  source.p_in=Mea.y[2];
  source.h_in=Mea.y[6];
  sink.p_in=Mea.y[1];

  connect(bulb.y,txv.p_b);
  connect(source.ports[1],txv.port_a);
  connect(txv.port_b,sink.ports[1]);

  annotation (experiment(
      StartTime=2000,
      StopTime=5000,
      Tolerance=0.001,
      __Dymola_Algorithm="Dassl"));
end TestTXV;
