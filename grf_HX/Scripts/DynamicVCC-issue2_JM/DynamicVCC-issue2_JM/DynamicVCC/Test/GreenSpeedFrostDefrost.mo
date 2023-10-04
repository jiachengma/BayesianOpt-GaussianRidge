within DynamicVCC.Test;
model GreenSpeedFrostDefrost "Greenspeed cycling of frosting and defrosting validaiton"
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

  parameter SI.CoefficientOfHeatTransfer alpha_nominal_ID=8000;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_OD=5000;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_SinglePhase=1000;
  parameter SI.LewisNumber Le=0.854 "Lewis number";
  parameter Integer N_pipe=5; //Discharge pipe line

/***************************** Initialization *******************************/

  parameter Real p_init_ID[Ncell]=fill(7.68e5,Ncell);
  parameter Real h_init_ID[Ncell]=fill(4.41e5,Ncell);
  parameter Real Tt_init_ID[Ncell]=fill(291.4,Ncell);
  parameter SI.MassFlowRate m_dot_air_ini_ID[Ncell]=fill(1e-3/Ncell,Ncell);
  parameter SI.Temperature T_a_init_ID[Ncell]=fill(292,Ncell);
  parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.004,Ncell);

  parameter Real p_init_OD[Ncell]=fill(7.6e5,Ncell);
  parameter Real h_init_OD[Ncell]=fill(4.22e5,Ncell);
  parameter Real Tt_init_OD[Ncell]=fill(272,Ncell);
  parameter SI.MassFlowRate m_dot_air_ini_OD[Ncell]=fill(1e-3/Ncell,Ncell);
  parameter SI.Temperature T_a_init_OD[Ncell]=fill(272,Ncell);
  parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.0023,Ncell);

  parameter SI.Frequency speed_init=1;
  parameter SI.Power Pwr_init=1;
  parameter Real opening_init=0.1;
  parameter SI.Pressure dp_init=23e5;
  parameter SI.MassFlowRate m_flows_init_1[Ncell+1]=fill(1e-3,Ncell+1);
  parameter SI.Volume V_f_init=0.0029;

  parameter SI.Thickness x_f_init=1e-5;
  parameter SI.Density rho_f_init=30;

/***************************** Numerical *******************************/
  import DynamicVCC.Components.Types.ModelStructure;
  import DynamicVCC.Components.Types.DifferentialState;

  parameter DifferentialState differentialState=DifferentialState.ph;
  parameter Boolean SteadyState_init=false;
  parameter Boolean EnableReverseFlow=true;
  parameter Boolean SteadyStateMomentum=false;
  parameter Boolean useLumpedPressure=false;
  parameter Boolean use_I_flows=true;
  parameter ModelStructure modelStructure=ModelStructure.av_vb;

/***************************** Heat Transfer and Pressure Drop *******************************/

  replaceable model Condensation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_ID);

  replaceable model Evaporation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_OD);

  replaceable model SinglePhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.SinglePhase_Gnielinski (
  alpha0=alpha_nominal_SinglePhase);

  replaceable model HeatTransfer_1_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
    redeclare model LiquidZone=SinglePhase,
    redeclare model VaporZone=SinglePhase,
    redeclare model TwoPhaseZone=Evaporation);

  replaceable model HeatTransfer_2_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_Wang (
    Nrow=Geo_OD.N_tuberow,
    pf=Geo_OD.pf,
    pl=Geo_OD.pl,
    pt=Geo_OD.pt,
    t_fin=Geo_OD.t_fin,
    final alpha0=45);

   replaceable model FreeConvection_OD=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_FreeConvection (
    final Const_alpha=false,
    final alpha0=1);

   replaceable model Friction_1_OD=DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (       f0=0.1,dp_nominal=1e3,m_flow_nominal=0.048);

   replaceable model Friction_2_OD = DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
     f0=0.12,m_flow_nominal=1)
     "Outdoor air side pressure drop";

   replaceable model HeatTransfer_1_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
      redeclare model LiquidZone = SinglePhase,
      redeclare model VaporZone = SinglePhase,
      redeclare model TwoPhaseZone = Condensation);

  replaceable model Friction_1_ID=DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (       f0=0.5,dp_nominal=1e2,m_flow_nominal=0.048);

 replaceable model HeatTransfer_2_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_Wang (
  Nrow=Geo_ID.N_tuberow,
  pf=Geo_ID.pf,
  pl=Geo_ID.pl,
  pt=Geo_ID.pt,
  t_fin=Geo_ID.t_fin,
  final alpha0=45)
  "Indoor air side heat transfer model";

 replaceable model FreeConvection_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_FreeConvection (
    final Const_alpha=false,
    final alpha0=1);

 replaceable model Friction_2_ID = DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
 f0=0.1,m_flow_nominal=1)
 "Indoor air side pressure drop";


/********************************* Control ******************************************/

  Modelica.Blocks.Continuous.FirstOrder Comp_speed(T=5,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=speed_init);
  Modelica.Blocks.Continuous.FirstOrder EXV_opening(T=1,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=opening_init);
  Modelica.Blocks.Continuous.FirstOrder RV_opening(T=3,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=1.0);
  Modelica.Blocks.Continuous.FirstOrder fan_OD_speed(T=10,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=1);
  Modelica.Blocks.Continuous.FirstOrder fan_ID_flow(T=1,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=sum(m_dot_air_ini_ID));

  Modelica.Blocks.Nonlinear.Limiter Comp_limiter(uMax=110,uMin=0.85);
  Modelica.Blocks.Nonlinear.Limiter exv_limiter(uMax=1.0,uMin=0.01);



  Modelica.Blocks.Sources.CombiTimeTable BC_ID(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_ID",fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/BC_ID.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable BC_OD(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_OD",fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/BC_OD.mat",columns=2:5);
  Modelica.Blocks.Sources.CombiTimeTable Mea_ID(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea_ID",fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/Mea_ID.mat",columns=2:6);
  Modelica.Blocks.Sources.CombiTimeTable Mea_OD(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea_OD",fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/Mea_OD.mat",columns=2:6);
  Modelica.Blocks.Sources.CombiTimeTable Comp(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Comp",fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/Comp.mat",columns=2:15);
  Modelica.Blocks.Sources.CombiTimeTable EXVopen(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="data",fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/EXVopen.mat",columns=2:2);

/***************************** Components *******************************/

  DynamicVCC.Components.Units.MassFlowDevices.Compressor.Efficiencies compressor(
    redeclare final package Medium = Medium_1,
    Vs=2.3765e-05,
    m_flow_nominal=0.047,
    speed_nominal=53,
    T_amb=T_amb,
    m_flow_init=m_flows_init_1[1],
    h_dis_init=h_init_ID[1],
    h_suc_init=h_init_OD[Ncell],
    p_dis_init=p_init_ID[1],
    p_suc_init=p_init_OD[Ncell],
    speed_init=speed_init,
    Pwr_init=Pwr_init);

  DynamicVCC.Components.Units.MassFlowDevices.Valve.ElectronicExpansionValve exv(
    redeclare package Medium = Medium_1,
    final Av=2.5447e-06,
    final EnableReverseFlow=EnableReverseFlow,
    opening_init=opening_init,
    dp_nominal=23e5,
    m_flow_nominal=0.047,
    dp_init=dp_init,
    m_flow_init=m_flows_init_1[1]);

  DynamicVCC.Components.Units.HX.FinTubeHX IndoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_ID,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_ID,
    redeclare final model FreeConvection = FreeConvection_ID,
    redeclare final model Friction_1 = Friction_1_ID,
    redeclare final model Friction_2 = Friction_2_ID,
    final differentialState=differentialState,
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
    w_a_init=w_a_init_ID,
    SteadyState_init=SteadyState_init,
    m_flows_a_init=m_dot_air_ini_ID,
    Ta_init=T_a_init_ID);

  DynamicVCC.Components.Units.HX.FrostDefrostHX.FinTubeHX_FrostDefrost OutdoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_OD,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_OD,
    redeclare final model Friction_1 = Friction_1_OD,
    redeclare final model Friction_2 = Friction_2_OD,
    final differentialState=differentialState,
    final modelStructure=modelStructure,
    final EnableReverseFlow=EnableReverseFlow,
    final useLumpedPressure=useLumpedPressure,
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
    w_a_init=w_a_init_OD,
    SteadyState_init=SteadyState_init,
    m_flows_a_init=m_dot_air_ini_OD,
    Ta_init=T_a_init_OD,
    final L_tube=Geo_OD.L_tube,
    final H_hx=Geo_OD.H_hx,
    final t_fin=Geo_OD.t_fin,
    final fin_meter=Geo_OD.fin_meter,
    final N_t_prow=Geo_OD.N_t_prow,
    final frostmode=frostmode,
    final T_amb=T_amb);

  DynamicVCC.Components.Units.MassFlowDevices.Valve.ReversingValve reversingValve(
    redeclare package Medium=Medium_1,
    m_flow_init=m_flows_init_1[1],
    p_dis_init=p_init_ID[1],
    p_suc_init=p_init_OD[Ncell],
    A_open=1.0e-4,
    dp_nominal=0.2e5,
    C_Hd=1,
    C_Hs=0.7,
    tau=1);

  DynamicVCC.Components.Units.MassFlowDevices.Accumulator accumulator(
    redeclare package Medium = Medium_1,
    m_flow_init=m_flows_init_1[1],
    V=0.0058,
    V_f_init=V_f_init,
    p_init=p_init_OD[Ncell],
    SteadyState_init=SteadyState_init);

  Modelica.Fluid.Sources.MassFlowSource_T fan_ID[Ncell](
    redeclare each final package Medium = Medium_2,
    each use_m_flow_in=true,
    each use_T_in=true,
    each use_X_in=true,
    each nPorts=1);

  DynamicVCC.Components.Units.MassFlowDevices.Fan fan_OD(
    final Ncell=Ncell,
    V_dot_init=1e-2);

  Modelica.Fluid.Sources.Boundary_pT sink_ID(redeclare package Medium = Medium_2, nPorts=Ncell);
  Modelica.Fluid.Sources.Boundary_pT sink_OD(redeclare package Medium = Medium_2, nPorts=Ncell);

  DynamicVCC.Components.Units.HX.Piping pipe(
    redeclare final package Medium_1 = Medium_1,
    redeclare final model HeatTransfer_1 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.ConstHTC (alpha0=100),
    redeclare final model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.HeatLossAmbient (alpha0=10),
    redeclare final model Friction_1 = Friction_1_OD,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final EnableReverseFlow=false,
    final SteadyState_init=SteadyState_init,
    final useLumpedPressure=useLumpedPressure,
    final use_I_flows=use_I_flows,
    final Ncell=N_pipe,
    final T_amb=T_amb,
    d_i=0.0110744,
    d_o=0.0127,
    L_tube=1.0922,
    p_init=fill(p_init_ID[1], N_pipe),
    h_init=fill(h_init_ID[1], N_pipe),
    Tt_init=fill(Tt_init_ID[1], N_pipe),
    m_flows_init=fill(m_flows_init_1[1], N_pipe + 1));

  Real frostmode;
  SI.Temperature T_amb;
  SI.Mass charge;
  SI.Energy energy "Energy conservation";


/***************************** Sensors ***********************************/

  DynamicVCC.Components.Units.Sensors.T_Superheat superheat(redeclare package Medium=Medium_1);
  DynamicVCC.Components.Units.Sensors.T_Superheat subcooling(redeclare package Medium=Medium_1);
  DynamicVCC.Components.Units.Sensors.T_Superheat superheat_IDVL(redeclare package Medium=Medium_1);
  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply_ID(
    redeclare package Medium=Medium_2,
    nPorts=Ncell);
  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply_OD(
    redeclare package Medium=Medium_2,
    nPorts=Ncell);

equation

  charge=IndoorCoil.charge+OutdoorCoil.charge+accumulator.charge+pipe.charge;
  energy=IndoorCoil.energy+OutdoorCoil.energy+accumulator.energy;

  fan_OD.T_in=BC_OD.y[1];
  fan_OD.X_in={BC_OD.y[3],1-BC_OD.y[3]};
  fan_OD_speed.u=if time<2457 then max(1,Comp.y[13]) elseif time<2805 then 1 elseif time<10185 then max(1,Comp.y[13]) elseif time<10825 then 1 else max(1,Comp.y[13]);

  connect(fan_OD_speed.y,fan_OD.speed);

  T_amb=BC_OD.y[1];

  fan_ID_flow.u=BC_ID.y[2];

  for i in 1:Ncell loop
     fan_ID[i].m_flow_in=fan_ID_flow.y/Ncell;
     fan_ID[i].T_in=BC_ID.y[1];
     fan_ID[i].X_in={BC_ID.y[3],1-BC_ID.y[3]};
  end for;

   Comp_speed.u=Comp.y[4];
   connect(Comp_speed.y,Comp_limiter.u);
   connect(compressor.speed,Comp_limiter.y);

   RV_opening.u=if time<2519 then 1.0 elseif time <2802 then 0.0 elseif time<10250 then 1.0 elseif time<10820 then 0.0 else 1.0;
   connect(RV_opening.y,reversingValve.opening_SP);

   EXV_opening.u=EXVopen.y[1];
   connect(EXV_opening.y,exv_limiter.u);
   connect(exv_limiter.y,exv.opening);

   frostmode=smooth(0,noEvent(if time<2440 then 1.0 elseif time<2441 then 1.0-(time-2440) elseif time<2810 then 0.0 elseif time<2811 then (time-2810) elseif time<10160 then 1.0 elseif time<10161 then 1.0-(time-10160) elseif time<10820 then 0.0 elseif time<10821 then time-10820 else 1.0));

  connect(superheat_IDVL.port, IndoorCoil.port_a1);
  connect(IndoorCoil.port_b1, exv.port_a);
  connect(exv.port_b, OutdoorCoil.port_a1);
  connect(fan_ID.ports[1], IndoorCoil.ports_a2);
  connect(fan_OD.ports, OutdoorCoil.ports_a2);
  connect(sink_ID.ports, IndoorCoil.ports_b2);
  connect(sink_OD.ports, OutdoorCoil.ports_b2);
  connect(superheat.port, OutdoorCoil.port_b1);
  connect(subcooling.port, IndoorCoil.port_b1);
  connect(T_supply_ID.ports, IndoorCoil.ports_b2);
  connect(T_supply_OD.ports, OutdoorCoil.ports_b2);
  connect(compressor.port_a, accumulator.port_b);
  connect(compressor.port_b, pipe.port_a1);
  connect(pipe.port_b1, reversingValve.port_a);
  connect(reversingValve.port_b, accumulator.port_a);
  connect(reversingValve.port_b2, IndoorCoil.port_a1);
  connect(OutdoorCoil.port_b1, reversingValve.port_a2);

  annotation (experiment(
      StartTime=615,
      StopTime=11830,
      Tolerance=0.01,
      __Dymola_Algorithm="Dassl"));
end GreenSpeedFrostDefrost;
