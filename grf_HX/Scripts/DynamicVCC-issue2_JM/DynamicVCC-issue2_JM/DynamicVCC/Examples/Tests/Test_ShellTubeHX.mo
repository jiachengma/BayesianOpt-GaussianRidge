within DynamicVCC.Examples.Tests;
model Test_ShellTubeHX "Test shell-and-tube HX model"
  extends Modelica.Icons.Example;

  //parameter Real u[5]={5e4,8e4,5e4,2e4,131146}; // condenser
  parameter Real u[5]={5e4,8e4,5e4,2e4,123886.44}; // evaporator
  output Real y[3];
  output Real y_mea[3];

  inner DynamicVCC.Components.System system(
  redeclare package Medium=Medium_1,
  T_max=330,
  T_min=250,
  m_flow_init=m_flow_init,
  m_flow_nominal=2.396,
  massDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
  energyDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
  momentumDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
  enableReverseFlow=true);

  package Medium_1=DynamicVCC.Media.R134a_NN;

  package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater (
  cp_const=4186.8,
  d_const=995);

  /* Initial conditions */

  // condenser
  /*
  parameter SI.Pressure p_init=9.2e5;
  parameter SI.SpecificEnthalpy h_init[Ncell]={425665,422357.093750000,420374.250000000,419181.812500000,418463.156250000,418029.500000000,417767.562500000,414636.218750000,409043.187500000,399053.312500000,381210.031250000,349339.531250000,292414.468750000,250096.875000000,247136.843750000};
  parameter SI.ThermodynamicTemperature Tt_init[Ncell]={309.448150634766,309.167724609375,309.000793457031,308.900909423828,308.840911865234,308.804748535156,308.782958984375,308.891448974609,308.885467529297,308.874816894531,308.855804443359,308.821807861328,308.761108398438,306.224029541016,303.171966552734};
  parameter SI.ThermodynamicTemperature Te_init[Ncell]={303.050842285156,304.492584228516,306.432006835938,307.517822265625,308.125732421875,308.466064453125,308.656616210938,308.763305664063,308.772247314453,308.787017822266,308.811492919922,308.852111816406,308.919677734375,309.032379150391,309.221282958984};
*/
  parameter SI.MassFlowRate m_flow_init=2.396;

  // evaporator
  parameter SI.Pressure p_init=3.91e5;
  parameter SI.SpecificEnthalpy h_init[Ncell]=linspace(2.51e5,4.06e5,Ncell);
  parameter SI.ThermodynamicTemperature Tt_init[Ncell]={281.610473632813,281.642669677734,281.686248779297,281.745239257813,281.825073242188,281.933166503906,282.079498291016,282.276245117188,284.541809082031,288.416748046875,288.960968017578,289.037475585938,289.048248291016,289.049743652344,289.049957275391};
  parameter SI.ThermodynamicTemperature Te_init[Ncell]={289.049987792969,289.049926757813,289.049407958984,289.045837402344,289.020477294922,288.840087890625,287.555786132813,285.978302001953,284.813385009766,283.952789306641,283.317047119141,282.847412109375,282.500457763672,282.244171142578,282.054809570313};


  /* Numerical */
  parameter Integer Ncell=15;

  // Numerics
  import DynamicVCC.Components.Types.ModelStructure;
  import DynamicVCC.Components.Types.DifferentialState;
  parameter ModelStructure modelStructure=ModelStructure.av_vb;
  parameter DifferentialState differentialState=DifferentialState.pdh;
  parameter Boolean useLumpedPressure=true;


  parameter SI.CoefficientOfHeatTransfer alpha_f=u[1];
  parameter SI.CoefficientOfHeatTransfer alpha_tp=u[2];
  parameter SI.CoefficientOfHeatTransfer alpha_g=u[3];
  parameter SI.CoefficientOfHeatTransfer alpha_w=u[4];
  parameter SI.HeatCapacity C_metalWall=u[5];

  /* Heat transfer */

  replaceable model HeatTransfer_1 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
  final alpha_f=alpha_f,
  final alpha_tp=alpha_tp,
  final alpha_g=alpha_g);

  replaceable model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
  alpha0=alpha_w);

  /* Friction */
  replaceable model FlowModel_1 = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
    final lambda0=0.01);

  replaceable model SlipRatio=DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.Homogeneous;

  //replaceable model RefFlow1D = DynamicVCC.Test.Chiller.RefFlow1D_UniformPressure;
  //replaceable model RefFlow1D = DynamicVCC.Components.Pipes.RefFlow1D;

  /* Condenser model */
/*
  DynamicVCC.Components.Units.HX.ShellTubeHX HX(
  redeclare final package Medium_1=Medium_1,
  redeclare final package Medium_2=Medium_2,
  redeclare final model HeatTransfer_1=HeatTransfer_1,
  redeclare final model HeatTransfer_2=HeatTransfer_2,
  redeclare final model SlipRatio=SlipRatio,
  redeclare final model FlowModel_1=FlowModel_1,
  final Ncell=Ncell,
  final C_tube=C_metalWall,
  As_1=23.9328,
  Ac_1=0.0652,
  L_1=2.4384,
  diameter_1=0.01905,
  As_2=19.5231,
  Ac_2=0.0311,
  L_2=2.4384,
  diameter_2=0.01554,
  M_metalWall=340.6385,
  cp_metalWall=385,
  modelStructure=modelStructure,
  differentialState=differentialState,
  useLumpedPressure=useLumpedPressure,
  p_a_start=p_init,
  h_init=h_init,
  Tt_init=Tt_init,
  Te_init=Te_init);
*/
  /* Evaporator model */

  DynamicVCC.Components.Units.HX.ShellTubeHX HX(
  redeclare final package Medium_1=Medium_1,
  redeclare final package Medium_2=Medium_2,
  redeclare final model HeatTransfer_1=HeatTransfer_1,
  redeclare final model HeatTransfer_2=HeatTransfer_2,
  redeclare final model SlipRatio=SlipRatio,
  redeclare final model FlowModel_1=FlowModel_1,
  final Ncell=Ncell,
  final C_tube=C_metalWall,
  As_1=22.3716,
  Ac_1=0.07627,
  L_1=2.4384,
  diameter_1=0.0196,
  As_2=18.331,
  Ac_2=0.03,
  L_2=2.4384,
  diameter_2=0.01606,
  M_metalWall=321.782954,
  cp_metalWall=385,
  modelStructure=modelStructure,
  differentialState=differentialState,
  useLumpedPressure=useLumpedPressure,
  p_a_start=p_init,
  h_init=h_init,
  Tt_init=Tt_init,
  Te_init=Te_init);


  /* Sources */
  Modelica.Fluid.Sources.MassFlowSource_h ref_source(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_h_in=true,
  use_m_flow_in=true);

  Modelica.Fluid.Sources.MassFlowSource_h ref_sink(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_m_flow_in=true);

  Modelica.Fluid.Sources.MassFlowSource_T water_source(
  redeclare package Medium=Medium_2,
  nPorts=1,
  use_T_in=true,
  use_m_flow_in=true);

  Modelica.Fluid.Sources.Boundary_pT water_sink(
  redeclare package Medium=Medium_2,
  nPorts=1,
  use_p_in=true);

  DynamicVCC.Components.Units.Sensors.T_Superheat superheat(
  redeclare final package Medium=Medium_1);


  Modelica.Blocks.Sources.CombiTimeTable BC_Cond(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Cond",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Cond.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable BC_Evap(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Evap",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Evap.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:9);

//Real y;
equation
  connect(ref_source.ports[1],HX.port_a1);
  connect(ref_sink.ports[1],HX.port_b1);
  connect(water_source.ports[1],HX.port_a2);
  connect(water_sink.ports[1],HX.port_b2);
  connect(HX.port_b1,superheat.port);
  water_sink.p_in=101325;

  // Outputs
  y[1]=HX.port_a1.p;
  y[2]=superheat.T;
  y[3]=HX.secondaryFluid.T[Ncell];

  y_mea[1]=Mea.y[2];
  y_mea[2]=Mea.y[8];
  y_mea[3]=Mea.y[4];


  // Boundary conditions
/*
  ref_source.m_flow_in=BC_Cond.y[1];
  ref_source.h_in=BC_Cond.y[2];
  ref_sink.m_flow_in=-BC_Evap.y[1];
  water_source.T_in=BC_Cond.y[3];
  water_source.m_flow_in=16.8;

  ref_source.m_flow_in=2.396393;
  ref_source.h_in=431209.618;
  ref_sink.m_flow_in=-2.396393;
  water_source.T_in=302.95;
  water_source.m_flow_in=16.8;


  ref_source.m_flow_in=BC_Evap.y[1];
  ref_source.h_in=BC_Evap.y[2];
  ref_sink.m_flow_in=-BC_Cond.y[1];
  water_source.T_in=BC_Evap.y[3];
  water_source.m_flow_in=BC_Evap.y[4];
  */
  ref_source.m_flow_in=2.396393;
  ref_source.h_in=2.43e5;
  ref_sink.m_flow_in=-2.396393;
  water_source.T_in=289.05;
  water_source.m_flow_in=13.7;


  annotation (experiment(
      StartTime=2000,
      StopTime=3500,
      __Dymola_Algorithm="Dassl"));
end Test_ShellTubeHX;
