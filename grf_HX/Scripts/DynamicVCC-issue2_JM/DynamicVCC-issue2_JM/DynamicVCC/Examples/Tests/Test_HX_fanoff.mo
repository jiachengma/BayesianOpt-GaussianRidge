within DynamicVCC.Examples.Tests;
model Test_HX_fanoff "HX under fan off condition"
 extends Modelica.Icons.Example;

   inner DynamicVCC.Components.System system(
  p_max=45e5,
  p_min=2e5,
  h_max=4.7e5,
  h_min=1.11e5,
  T_max=340,
  T_min=250,
  EnableReverseFlow=EnableReverseFlow);

  parameter Integer Ncell=30;

  //replaceable package Medium_CP=DynamicVCC.Media.CoolProp.R410a;
  replaceable package Medium_NN=DynamicVCC.Media.R410A_ph;
  package Medium_1=Medium_NN;
  package Medium_2=Modelica.Media.Air.MoistAir;

/***************************** HX Geometry *******************************/

  DynamicVCC.Examples.Tests.FinTubeGeo Geo_OD(
    D_o=0.0074,
    D_i=0.0068,
    cp_tube=385,
    rho_tube=8900,
    cp_fin=900,
    rho_fin=2700,
    L_tube=2.7026,
    N_tuberow=2,
    N_t_prow=44,
    N_circuits=8,
    pf=0.0012,
    pt=0.0216,
    pl=0.0187,
    fin_meter=787.4016,
    t_fin=9.9060e-05);

  DynamicVCC.Examples.Tests.FinTubeGeo Geo_ID(
    D_o=0.01,
    D_i=0.0086,
    cp_tube=900,
    rho_tube=2700,
    cp_fin=900,
    rho_fin=2700,
    L_tube=0.4521,
    N_tuberow=3,
    N_t_prow=28,
    N_circuits=6,
    pf=0.0016,
    pt=0.0254,
    pl=0.0191,
    fin_meter=570.8661,
    t_fin=1.0668e-04);

/***************************** Constants *******************************/

  parameter SI.AbsolutePressure p_scale=1e6;
  parameter SI.SpecificEnthalpy h_scale=1e5;
  parameter SI.MassFlowRate m_scale=1e-2;
  parameter SI.Temperature T_scale=273;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_ID=8000;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_OD=5000;
  parameter SI.LewisNumber Le=0.854 "Lewis number";

/***************************** Initialization *******************************/
// A steady state

      parameter Real p_init_ID[Ncell]=linspace(32.3e5,32.2e5,Ncell);
      parameter Real h_init_ID[Ncell]=linspace(4.34e5,2.5e5,Ncell);
      parameter Real Tt_init_ID[Ncell]=linspace(318,302,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_ID=0.435;
      parameter SI.Temperature T_a_init_ID[Ncell]=linspace(317,301,Ncell);
      parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.002,Ncell);

      parameter Real p_init_OD[Ncell]=linspace(12.6e5,11.9e5,Ncell);
      parameter Real h_init_OD[Ncell]=linspace(2.58e5,4.25e5,Ncell);
      parameter Real Tt_init_OD[Ncell]=linspace(288.7,287.7,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_OD=1.6;
      parameter SI.Temperature T_a_init_OD[Ncell]=linspace(288,287,Ncell);
      parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.002,Ncell);
      parameter SI.MassFlowRate m_flows_init_1[Ncell+1]=fill(0.048,Ncell+1);

// start-up
/*
  parameter Real p_init_ID[Ncell]=linspace(14.18,14.17,Ncell);
  parameter Real h_init_ID[Ncell]=fill(4.41e5,Ncell);
  parameter Real Tt_init_ID[Ncell]=fill(291.4,Ncell);
  parameter SI.MassFlowRate m_dot_air_ini_ID[Ncell]=fill(1e-3/Ncell,Ncell);
  parameter SI.Temperature T_a_init_ID[Ncell]=fill(292,Ncell);
  parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.004,Ncell);

  parameter Real p_init_OD[Ncell]=fill(14.18e5,Ncell);
  parameter Real h_init_OD[Ncell]=fill(4.21e5,Ncell);
  parameter Real Tt_init_OD[Ncell]=fill(293,Ncell);
  parameter SI.MassFlowRate m_dot_air_ini_OD=1e-3;
  parameter SI.Temperature T_a_init_OD[Ncell]=fill(293,Ncell);
  parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.0022,Ncell);
  parameter SI.MassFlowRate m_flows_init_1[Ncell+1]=fill(1e-5,Ncell+1);
*/
/***************************** Numerical *******************************/
  import DynamicVCC.Components.Types.ModelStructure;
  parameter Boolean SteadyState_init=false;
  parameter Boolean EnableReverseFlow=true;
  parameter Boolean SteadyStateMomentum=false;
  parameter Boolean useLumpedPressure=false;
  parameter Boolean use_I_flows=true;
  parameter ModelStructure modelStructure=ModelStructure.av_vb;

/***************************** Heat Transfer and Pressure Drop *******************************/

/*
  replaceable model Condensation =DynamicVCC.Components.Pipes_newcell.BaseClasses.HeatTransfer.Correlations.Condensation_Shah (
  pc=4901200,
  alpha0=alpha_nominal_ID);
*/

  replaceable model Condensation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_ID);

  replaceable model Evaporation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_OD);

  replaceable model Singlephase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.SinglePhase_Gnielinski (
  alpha0=alpha_nominal_OD);

 /*
  replaceable model Singlephase =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  alpha0=500);
 */
  replaceable model HeatTransfer_1_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=Singlephase,
  redeclare model VaporZone=Singlephase,
  redeclare model TwoPhaseZone=Evaporation);

 replaceable model HeatTransfer_2_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_Wang (
  Nrow=Geo_OD.N_tuberow,
  pf=Geo_OD.pf,
  pl=Geo_OD.pl,
  pt=Geo_OD.pt,
  t_fin=Geo_OD.t_fin,
  final alpha0=45);

  replaceable model FreeConvection_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_FreeConvection (
    final Const_alpha=true,
    final alpha0=0.1);

/*
  replaceable model HeatTransfer_2_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_ConstCoefficient (
    final alpha0=50);
 */

 replaceable model Friction_1_OD=DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (       f0=0.1);

 replaceable model Friction_2_OD = DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
 f0=0.12)
 "Outdoor air side pressure drop";

  replaceable model HeatTransfer_1_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=Singlephase,
  redeclare model VaporZone=Singlephase,
  redeclare model TwoPhaseZone=Condensation);

 //replaceable model Friction_1_ID=DynamicVCC.Components.Pipes_newcell.BaseClasses.Friction.Correlations.TwoPhase_SinglePhase;
 replaceable model Friction_1_ID=DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (       f0=0.5);

 replaceable model HeatTransfer_2_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_Wang (
  Nrow=Geo_ID.N_tuberow,
  pf=Geo_ID.pf,
  pl=Geo_ID.pl,
  pt=Geo_ID.pt,
  t_fin=Geo_ID.t_fin,
  final alpha0=45)
  "Indoor air side heat transfer model";

/*
  replaceable model Friction_2_ID=DynamicVCC.Components.Pipes_newcell.BaseClasses.Friction.AirCoilDP_Wang (
      Nt=Geo_ID.N_tuberow,
      D=Geo_ID.D_o,
      pf=Geo_ID.pf,
      pl=Geo_ID.pl,
      pt=Geo_ID.pt,
      t_fin=Geo_ID.t_fin,
      final dp_nominal=20);
*/

 replaceable model Friction_2_ID = DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
 f0=0.1)
 "Indoor air side pressure drop";

 replaceable model FreeConvection =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_FreeConvection (
 final alpha0=0.1,
 final Const_alpha=true);

/***************************** Components *******************************/

/*
  DynamicVCC.Components.Units.HX.FinTubeHX IndoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_ID,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_ID,
    redeclare final model FreeConvection = FreeConvection,
    redeclare final model Friction_1 = Friction_1_ID,
    redeclare final model Friction_2 = Friction_2_ID,
    final modelStructure=modelStructure,
    final EnableReverseFlow=EnableReverseFlow,
    final SteadyStateMomentum=SteadyStateMomentum,
    final useLumpedPressure=useLumpedPressure,
    final use_I_flows=use_I_flows,
    final Ncell=Ncell,
    As_1=Geo_ID.HTA_r,
    Ac_1=Geo_ID.Ac_r,
    L_1=Geo_ID.L_circuit,
    M_metalWall=Geo_ID.M_tube,
    cp_metalWall=Geo_ID.cp_tube,
    diameter_2=Geo_ID.D_o,
    diameter_1=Geo_ID.D_i,
    L_2=Geo_ID.pl,
    Eta_fin_overall=Geo_ID.Eta_fin_overall,
    Le=Le,
    As_2=Geo_ID.HTA_e,
    Ac_2=Geo_ID.Ac_e,
    M_fin=Geo_ID.M_fin,
    cp_fin=Geo_ID.cp_fin,
    p_init=p_init_ID,
    h_init=h_init_ID,
    Tt_init=Tt_init_ID,
    m_flows_init=m_flows_init_1,
    p_scale=p_scale,
    h_scale=h_scale,
    T_scale=T_scale,
    w_a_init=w_a_init_ID,
    SteadyState_init=SteadyState_init,
    m_flows_a_init=fill(m_dot_air_ini_ID/Ncell,Ncell),
    Ta_init=T_a_init_ID) annotation (Placement(transformation(extent={{16,30},{-56,102}})));

*/

  Modelica.Fluid.Sources.MassFlowSource_T fan[Ncell](
  redeclare each final package Medium=Medium_2,
  each use_m_flow_in=true,
  each use_T_in=true,
  each use_X_in=true,
  each nPorts=1);

  Modelica.Fluid.Sources.Boundary_pT sink_air(
  redeclare package Medium=Medium_2,
  nPorts=Ncell);

  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply(redeclare package Medium=Medium_2,
  nPorts=Ncell);

  DynamicVCC.Components.Units.Sensors.T_Superheat sensor_sh(redeclare package Medium=Medium_1);

  Components.Units.HX.FinTubeHX HX(
  redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_OD,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_OD,
    redeclare final model Friction_1 = Friction_1_OD,
    redeclare final model Friction_2 = Friction_2_OD,
    final modelStructure=modelStructure,
    final EnableReverseFlow=EnableReverseFlow,
    final useLumpedPressure=useLumpedPressure,
    final SteadyStateMomentum=SteadyStateMomentum,
    final use_I_flows=use_I_flows,
    final Ncell=Ncell,
    As_1=Geo_OD.HTA_r,
    Ac_1=Geo_OD.Ac_r,
    L_1=Geo_OD.L_circuit,
    M_metalWall=Geo_OD.M_tube,
    cp_metalWall=Geo_OD.cp_tube,
    diameter_2=Geo_OD.D_o,
    diameter_1=Geo_OD.D_i,
    L_2=Geo_OD.pl,
    Eta_fin_overall=Geo_OD.Eta_fin_overall,
    Le=Le,
    As_2=Geo_OD.HTA_e,
    Ac_2=Geo_OD.Ac_e,
    M_fin=Geo_OD.M_fin,
    cp_fin=Geo_OD.cp_fin,
    p_init=p_init_OD,
    h_init=h_init_OD,
    Tt_init=Tt_init_OD,
    m_flows_init=m_flows_init_1,
    p_scale=p_scale,
    h_scale=h_scale,
    m_scale=m_scale,
    T_scale=T_scale,
    w_a_init=w_a_init_OD,
    SteadyState_init=SteadyState_init,
    m_flows_a_init=fill(m_dot_air_ini_OD/Ncell,Ncell),
    Ta_init=T_a_init_OD)  annotation (Placement(transformation(extent={{-12,4},{38,54}})));

  Modelica.Fluid.Sources.MassFlowSource_h source(
    redeclare final package Medium=Medium_1,
    use_m_flow_in=true,
    use_h_in=true,
    nPorts=1) annotation (Placement(transformation(extent={{-64,16},{-44,36}})));
  Modelica.Fluid.Sources.MassFlowSource_h sink(
  redeclare final package Medium=Medium_1,
  use_m_flow_in=true, nPorts=1) annotation (Placement(transformation(extent={{76,18},{56,38}})));

equation
  /*
    source.m_flow_in=1e-5;
    source.h_in=4.21e5;
    sink.m_flow_in=-1e-5;
   */
    source.m_flow_in=0.048;
    source.h_in=2.54e5;
    sink.m_flow_in=-0.048;

  /*
   for i in 1:Ncell loop
     fan[i].m_flow_in=m_OD_air/Ncell;
     fan[i].X_in={BC_OD.y[3],1-BC_OD.y[3]};
     fan[i].T_in=BC_OD.y[1];
   end for;
   */

  for i in 1:Ncell loop
     fan[i].m_flow_in=1.6/Ncell;
     fan[i].X_in={0.0022,1-0.0022};
     fan[i].T_in=293;
  end for;

/*
  for i in 1:Ncell loop
   fan[i].m_flow_in=1e-3/Ncell;
   fan[i].T_in=293;
   fan[i].X_in={0.0022,1-0.0022};
  end for;
 */
  connect(fan.ports[1],HX.ports_a2);
  connect(HX.ports_b2,sink_air.ports);
  connect(HX.ports_b2,T_supply.ports);
  connect(HX.port_b1,sensor_sh.port);

  connect(HX.port_b1, sink.ports[1]) annotation (Line(points={{38,29},{38,28},{56,28}}, color={0,127,255}));
  connect(source.ports[1], HX.port_a1) annotation (Line(points={{-44,26},{-18,26},{-18,29},{-12,29}},     color={0,127,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=2000,
      Tolerance=0.001,
      __Dymola_Algorithm="Dassl"));
end Test_HX_fanoff;
