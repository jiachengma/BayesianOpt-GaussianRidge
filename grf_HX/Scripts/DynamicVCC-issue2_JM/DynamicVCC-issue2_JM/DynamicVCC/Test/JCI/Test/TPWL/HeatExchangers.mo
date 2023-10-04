within DynamicVCC.Test.JCI.Test.TPWL;
partial model HeatExchangers

  DynamicVCC.Components.Units.HX.BaseClasses.FinTubeCoil Geo_OD(
  d_o=0.007,
  d_i=0.00614,
  cp_tube=385,
  rho_tube=8960,
  cp_fin=900,
  rho_fin=2710,
  L_tube=2.5104344,
  N_row=2,
  N_prow=44,
  N_circuits=8,
  pf=0.00131205,
  pl=0.015875,
  pt=0.020320,
  t_fin=0.00009906,
  Eta_fin_overall=1,
  Ac_e=1.349,
  HTA_e=87.7181,
  M_fin=11.774,
  HTA_r=4.261375,
  Ac_r=0.00023687,
  M_tube=17.56803779);







  DynamicVCC.Components.Units.HX.BaseClasses.FinTubeCoil Geo_ID(
  d_o=0.00955,
  d_i=0.007518,
  cp_tube=900,
  rho_tube=2700,
  cp_fin=893,
  rho_fin=2720,
  L_tube=0.470154,
  N_row=3,
  N_prow=28,
  N_circuits=8,
  pf=0.0018398,
  pl=0.021996,
  pt=0.0254,
  t_fin=0.000114,
  Eta_fin_overall=1,
  Ac_e=0.18151322,
  HTA_e=18.19055,
  M_fin=2.820263986,
  HTA_r=0.9327636549,
  Ac_r=0.000355127,
  M_tube=2.904553663);

  replaceable package Medium_1=DynamicVCC.Media.R410a_NN;
  replaceable package Medium_2=Modelica.Media.Air.MoistAir;

  inner DynamicVCC.Components.System system(
  redeclare package Medium=Medium_1,
  T_max=330,
  T_min=250,
  m_flow_init=m_flow_init,
  m_flow_nominal=0.05,
  massDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
  energyDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
  momentumDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
  enableReverseFlow=true);

  parameter Integer Ncell=20;

  // Initial conditions
  parameter Medium_1.AbsolutePressure p_a_start_OD=23.7e5;
  parameter Medium_1.AbsolutePressure p_b_start_OD=p_a_start_OD-5.8343e4;
  parameter Medium_1.SpecificEnthalpy h_init_OD[Ncell]=linspace(4.7e5,2.55e5,Ncell);
  parameter SI.Temperature Tt_init_OD[Ncell]=fill(310,Ncell);

  parameter Medium_1.AbsolutePressure p_a_start_ID=p_b_start_ID+7.4601e4;
  parameter Medium_1.AbsolutePressure p_b_start_ID=9.3e5;
  parameter Medium_1.SpecificEnthalpy h_init_ID[Ncell]=linspace(2.55e5,4.34e5,Ncell);
  parameter SI.Temperature Tt_init_ID[Ncell]=fill(278,Ncell);
  parameter Medium_1.MassFlowRate m_flow_init=0.05;

  // Numerics
  import DynamicVCC.Components.Types.ModelStructure;
  import DynamicVCC.Components.Types.DifferentialState;
  parameter ModelStructure modelStructure=ModelStructure.av_vb;
  parameter DifferentialState differentialState=DifferentialState.pdh;
  parameter Boolean useLumpedPressure=false;

  //parameter Real u_OD[5]={1e3,1e3,8e3,51,17360.3};
  parameter Real u_OD[5]={750,750,6000,38.25,13020.225};
  parameter Real u_ID[4]={1e3,2e4,80,5132.6};
  parameter Real u_ss[:]=cat(1,u_OD,u_ID);
  //parameter Real u[9]={1e3,1e3,1e4,150,17360.3,500,5e3,70,5132.6};

  parameter SI.CoefficientOfHeatTransfer alpha_f_OD=u_ss[1];
  parameter SI.CoefficientOfHeatTransfer alpha_g_OD=u_ss[2];
  parameter SI.CoefficientOfHeatTransfer alpha_tp_OD=u_ss[3];
  parameter SI.CoefficientOfHeatTransfer alpha_a_OD=u_ss[4];
  parameter SI.HeatCapacity C_MetalWall_OD=u_ss[5];
/*
  parameter SI.CoefficientOfHeatTransfer alpha_f_OD=u_OD[1];
  parameter SI.CoefficientOfHeatTransfer alpha_g_OD=u_OD[2];
  parameter SI.CoefficientOfHeatTransfer alpha_tp_OD=u_OD[3];
  parameter SI.CoefficientOfHeatTransfer alpha_a_OD=u_OD[4];
  parameter SI.HeatCapacity C_MetalWall_OD=u_OD[5];
*/

  parameter SI.CoefficientOfHeatTransfer alpha_f_ID=1000;
  parameter SI.CoefficientOfHeatTransfer alpha_g_ID=u_ss[6];
  parameter SI.CoefficientOfHeatTransfer alpha_tp_ID=u_ss[7];
  parameter SI.CoefficientOfHeatTransfer alpha_a_ID=u_ss[8];
  parameter SI.HeatCapacity C_MetalWall_ID=u_ss[9];

 /*
  parameter SI.CoefficientOfHeatTransfer alpha_f_ID=1e3;
  parameter SI.CoefficientOfHeatTransfer alpha_g_ID=u_ID[1];
  parameter SI.CoefficientOfHeatTransfer alpha_tp_ID=u_ID[2];
  parameter SI.CoefficientOfHeatTransfer alpha_a_ID=u_ID[3];
  parameter SI.HeatCapacity C_MetalWall_ID=u_ID[4];
  */

  parameter SI.CoefficientOfFriction lambda_OD=0.02;
  parameter SI.CoefficientOfFriction lambda_ID=0.02;

  parameter SI.LewisNumber Le=0.854 "Lewis number";

  // Heat transfer and pressure drop

  replaceable model HeatTransfer_1_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
  final alpha_f=alpha_f_OD,
  final alpha_tp=alpha_tp_OD,
  final alpha_g=alpha_g_OD);

 /*
  replaceable model HeatTransfer_1_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.CorrelationPhaseChange (
  redeclare final model TwoPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
  final alpha_nominal=alpha_tp_OD,final m_flow_nominal=system.m_flow_nominal,final k=1, final b=0.5, final crossAreas=fill(Geo_OD.Ac_r,Ncell)),
  redeclare final model VaporPhase =
  DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
  final alpha_nominal=alpha_g_OD,final m_flow_nominal=system.m_flow_nominal,final k=1, final b=0.5,final crossAreas=fill(Geo_OD.Ac_r,Ncell)),
  redeclare final model LiquidPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
  final alpha_nominal=alpha_f_OD,final m_flow_nominal=system.m_flow_nominal,final k=1, final b=0.5,final crossAreas=fill(Geo_OD.Ac_r,Ncell)));
*/
  replaceable model HeatTransfer_2_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
  final alpha0=alpha_a_OD);

/*
  replaceable model HeatTransfer_2_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.PlainFinAndTube_Wang (
  nRow=Geo_OD.N_row,
  pf=Geo_OD.pf,
  pt=Geo_OD.pt,
  pl=Geo_OD.pl,
  t_fin=Geo_OD.t_fin,
  Dh=4*Geo_OD.Ac_e*Geo_OD.N_row*Geo_OD.pl/Geo_OD.HTA_e,
  alpha0=50);
  */

  replaceable model FreeConvection_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
  final alpha0=0.1);

  replaceable model FlowModel_1_OD = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=lambda_OD);

/*
  replaceable model FlowModel_1_OD=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
  final dp_nominal=5.8343e4/Ncell,
  final k=0.9808,
  final b=0.7512);
  */



  replaceable model Friction_2_OD=DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.ConstantFriction (
  f0=0.1);

 /*
  replaceable model Friction_2_OD =DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.PlainFinAndTube_Wang (
  nRow=Geo_OD.N_row,
  pf=Geo_OD.pf,
  pt=Geo_OD.pt,
  pl=Geo_OD.pl,
  t_fin=Geo_OD.t_fin,
  Dh=4*Geo_OD.Ac_e*Geo_OD.N_row*Geo_OD.pl/Geo_OD.HTA_e,
  f0=0.1);
*/

  replaceable model HeatTransfer_1_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
      final alpha_f=alpha_f_ID,
      final alpha_tp=alpha_tp_ID,
      final alpha_g=alpha_g_ID);

/*
  replaceable model HeatTransfer_1_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.CorrelationPhaseChange (
    redeclare final model TwoPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
    final alpha_nominal=alpha_tp_ID,final m_flow_nominal=system.m_flow_nominal,final k=1, final b=0.5,final crossAreas=fill(Geo_ID.Ac_r,Ncell)),
    redeclare final model VaporPhase =
    DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
    final alpha_nominal=alpha_g_ID,final m_flow_nominal=system.m_flow_nominal,final k=1, final b=0.5,final crossAreas=fill(Geo_ID.Ac_r,Ncell)),
    redeclare final model LiquidPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
    final alpha_nominal=alpha_f_ID,final m_flow_nominal=system.m_flow_nominal,final k=1, final b=0.5,final crossAreas=fill(Geo_ID.Ac_r,Ncell)));
    */

  replaceable model HeatTransfer_2_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
  final alpha0=alpha_a_ID);

 /*
  replaceable model HeatTransfer_2_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.PlainFinAndTube_Wang (
  nRow=Geo_ID.N_row,
  pf=Geo_ID.pf,
  pt=Geo_ID.pt,
  pl=Geo_ID.pl,
  t_fin=Geo_ID.t_fin,
  Dh=4*Geo_ID.Ac_e*Geo_ID.N_row*Geo_ID.pl/Geo_ID.HTA_e,
  alpha0=50);
 */
  replaceable model FreeConvection_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
  final alpha0=0.1);


  replaceable model FlowModel_1_ID=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=lambda_ID);

 /*
  replaceable model FlowModel_1_ID =DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
  final dp_nominal=7.4601e4/Ncell,
  final k=1.0065,
  final b=2.1559);
 */
  replaceable model Friction_2_ID =
  DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.ConstantFriction (
  f0=0.1);

  replaceable model SlipRatio=DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.Homogeneous;


  DynamicVCC.Components.Units.HX.FinTubeHX OutdoorCoil(
  redeclare final package Medium_1=Medium_1,
  redeclare final package Medium_2=Medium_2,
  redeclare final model HeatTransfer_1=HeatTransfer_1_OD,
  redeclare final model HeatTransfer_2=HeatTransfer_2_OD,
  redeclare final model FreeConvection=FreeConvection_OD,
  redeclare final model SlipRatio=SlipRatio,
  redeclare final model FlowModel_1=FlowModel_1_OD,
  redeclare final model FlowModel_2=Friction_2_OD,
  final modelStructure=modelStructure,
  final differentialState=differentialState,
  final useLumpedPressure=useLumpedPressure,
  final Ncell=Ncell,
  final C_FinTube=C_MetalWall_OD,
  As_1=Geo_OD.HTA_r,
  Ac_1=0.000355310359209821,
  L_1=18.4098522666667,
  diameter_1=Geo_OD.d_i,
  Ac_2=Geo_OD.Ac_e,
  As_2=Geo_OD.HTA_e,
  diameter_2=Geo_OD.d_o,
  L_2=Geo_OD.L_fin,
  p_a_start=p_a_start_OD,
  p_b_start=p_b_start_OD,
  h_init=h_init_OD,
  Tt_init=Tt_init_OD) annotation (Placement(transformation(extent={{-22,22},{22,-22}},
        rotation=180,
        origin={0,60})));


  DynamicVCC.Components.Units.HX.FinTubeHX IndoorCoil(
  redeclare final package Medium_1=Medium_1,
  redeclare final package Medium_2=Medium_2,
  redeclare final model HeatTransfer_1=HeatTransfer_1_ID,
  redeclare final model HeatTransfer_2=HeatTransfer_2_ID,
  redeclare final model FreeConvection=FreeConvection_ID,
  redeclare final model SlipRatio=SlipRatio,
  redeclare final model FlowModel_1=FlowModel_1_ID,
  redeclare final model FlowModel_2=Friction_2_ID,
  final modelStructure=modelStructure,
  final differentialState=differentialState,
  final useLumpedPressure=useLumpedPressure,
  final Ncell=Ncell,
  final C_FinTube=C_MetalWall_ID,
  As_1=Geo_ID.HTA_r,
  Ac_1=Geo_ID.Ac_r,
  L_1=Geo_ID.L_circuit,
  diameter_1=Geo_ID.d_i,
  Ac_2=Geo_ID.Ac_e,
  As_2=Geo_ID.HTA_e,
  diameter_2=Geo_ID.d_o,
  L_2=Geo_ID.L_fin,
  final p_a_start=p_a_start_ID,
  final p_b_start=p_b_start_ID,
  h_init=h_init_ID,
  Tt_init=Tt_init_ID) annotation (Placement(transformation(extent={{-22,-82},{22,-38}})));



  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=1,
      StopTime=1000,
      __Dymola_Algorithm="Dassl"));
end HeatExchangers;
