within DynamicVCC.Examples;
model AirSourceHeatPump "An R410A air-source heat pump system model"
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
  parameter SI.Temperature T_scale=273;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_ID=8000;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_OD=5000;
  parameter SI.LewisNumber Le=0.854 "Lewis number";
  parameter Integer N_pipe=5; //Discharge pipe line
  parameter SI.Time samplePeriod_exv=90;
  parameter SI.Time startTime=0 "Simulation start time";

/***************************** Initialization *******************************/

      parameter Real p_init_ID[Ncell]=linspace(33.6e5,33.5e5,Ncell);
      parameter Real h_init_ID[Ncell]=linspace(4.35e5,2.5e5,Ncell);
      parameter Real Tt_init_ID[Ncell]=linspace(318,302,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_ID[Ncell]=fill(0.435/Ncell,Ncell);
      parameter SI.Temperature T_a_init_ID[Ncell]=linspace(317,301,Ncell);
      parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.002,Ncell);

      parameter Real p_init_OD[Ncell]=linspace(12.7e5,11.7e5,Ncell);
      parameter Real h_init_OD[Ncell]=linspace(2.58e5,4.25e5,Ncell);
      parameter Real Tt_init_OD[Ncell]=linspace(288.7,287.7,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_OD[Ncell]=fill(1.6/Ncell,Ncell);
      parameter SI.Temperature T_a_init_OD[Ncell]=linspace(288,287,Ncell);
      parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.002,Ncell);

      parameter SI.Frequency speed_init=54.5;
      parameter SI.Power Pwr_init=1993;
      parameter Real opening_init=0.16;
      parameter SI.Pressure dp_init=20.9e5;
      parameter SI.MassFlowRate m_flows_init_1[Ncell+1]=fill(0.048,Ncell+1);
      parameter SI.Volume V_f_init=0.0028;

/***************************** Numerical *******************************/
  import DynamicVCC.Components.Types.ModelStructure;
  parameter Boolean SteadyState_init=false;
  parameter Boolean EnableReverseFlow=true;
  parameter Boolean useLumpedPressure=false;
  parameter Boolean use_I_flows=true;
  parameter DynamicVCC.Components.Types.ModelStructure modelStructure=DynamicVCC.Components.Types.ModelStructure.av_vb;

/***************************** Heat Transfer and Pressure Drop *******************************/

/*
  replaceable model Condensation =DynamicVCC.Components.Pipes_newcell.BaseClasses.HeatTransfer.Correlations.Condensation_Shah (
  pc=4901200,
  alpha0=alpha_nominal_ID);
*/

  replaceable model Condensation =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_ID);

  replaceable model Evaporation =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_OD);

  replaceable model Singlephase =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.SinglePhase_Gnielinski (
  alpha0=alpha_nominal_OD);


  replaceable model HeatTransfer_1_OD =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=Singlephase,
  redeclare model VaporZone=Singlephase,
  redeclare model TwoPhaseZone=Evaporation);

 replaceable model HeatTransfer_2_OD =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_Wang (
  Nrow=Geo_OD.N_tuberow,
  pf=Geo_OD.pf,
  pl=Geo_OD.pl,
  pt=Geo_OD.pt,
  t_fin=Geo_OD.t_fin,
  final alpha0=45);

  replaceable model FreeConvection_OD =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_FreeConvection (
    final Const_alpha=true,
    final alpha0=0.1);


 replaceable model Friction_1_OD =
      DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (                                  f0=0.1);

 replaceable model Friction_2_OD =
      DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
 f0=0.12)
 "Outdoor air side pressure drop";

  replaceable model HeatTransfer_1_ID =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=Singlephase,
  redeclare model VaporZone=Singlephase,
  redeclare model TwoPhaseZone=Condensation);

  replaceable model Friction_1_ID =
      DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (                                   f0=0.5);

 replaceable model HeatTransfer_2_ID =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_Wang (
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

 replaceable model Friction_2_ID =
      DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
 f0=0.1)
 "Indoor air side pressure drop";

/********************************* Control ******************************************/

  Modelica.Blocks.Continuous.FirstOrder fan_ID_flow(T=1,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=sum(m_dot_air_ini_ID));

  /*
// EXV controller
  Modelica.Blocks.Discrete.Sampler sampler_sh(samplePeriod=samplePeriod_exv,startTime=startTime);
  Modelica.Blocks.Discrete.ZeroOrderHold ZOH_exv(samplePeriod=samplePeriod_exv,startTime=startTime);
  Modelica.Blocks.Continuous.PI exv_controller(
  k=-2e-3,
  T=100,
  y_start=opening_init,
  x_start=-80,
  initType=Modelica.Blocks.Types.Init.InitialOutput);
  */

/***************************** Components *******************************/

  DynamicVCC.Components.Units.MassFlowDevices.Compressor.Efficiencies compressor(
    redeclare final package Medium = Medium_1,
    Vs=2.3765e-05,
    T_amb=T_amb,
    m_flow_init=m_flows_init_1[1],
    h_dis_init=h_init_ID[1],
    h_suc_init=h_init_OD[Ncell],
    p_dis_init=p_init_ID[1],
    p_suc_init=p_init_OD[Ncell],
    speed_init=speed_init,
    Pwr_init=Pwr_init) annotation (Placement(transformation(
        extent={{-26,-26},{26,26}},
        rotation=90,
        origin={148,18})));

  DynamicVCC.Components.Units.MassFlowDevices.Valve.ElectronicExpansionValve exv(
    redeclare package Medium = Medium_1,
    final Av=2.5447e-06,
    final EnableReverseFlow=EnableReverseFlow,
    p_a_init=p_init_ID[Ncell],
    p_b_init=p_init_OD[1],
    opening_init=opening_init,
    dp_nominal=23e5,
    dp_init=dp_init,
    m_flow_init=m_flows_init_1[1]) annotation (Placement(transformation(
        extent={{14,14},{-14,-14}},
        rotation=90,
        origin={-134,22})));

  DynamicVCC.Components.Units.HX.FinTubeHX IndoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_ID,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_ID,
    redeclare final model Friction_1 = Friction_1_ID,
    redeclare final model Friction_2 = Friction_2_ID,
    final modelStructure=modelStructure,
    final EnableReverseFlow=EnableReverseFlow,
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
    m_flows_a_init=m_dot_air_ini_ID,
    Ta_init=T_a_init_ID) annotation (Placement(transformation(extent={{16,30},{-56,102}})));

    DynamicVCC.Components.Units.HX.FinTubeHX OutdoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_OD,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_OD,
    redeclare final model Friction_1 = Friction_1_OD,
    redeclare final model Friction_2 = Friction_2_OD,
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
    p_scale=p_scale,
    h_scale=h_scale,
    T_scale=T_scale,
    w_a_init=w_a_init_OD,
    SteadyState_init=SteadyState_init,
    m_flows_a_init=m_dot_air_ini_OD,
    Ta_init=T_a_init_OD) annotation (Placement(transformation(extent={{-52,-92},{20,-20}})));

  DynamicVCC.Components.Units.MassFlowDevices.Valve.ReversingValve reversingValve(
    redeclare package Medium = Medium_1,
    m_flow_init=m_flows_init_1[1],
    p_dis_init=p_init_ID[1],
    p_suc_init=p_init_OD[Ncell],
    A_open=1.0e-4,
    dp_nominal=0.2e5,
    C_Hd=1,
    C_Hs=0.7,
    tau=1) annotation (Placement(transformation(
        extent={{-19,-19},{19,19}},
        rotation=-90,
        origin={79,-7})));

  DynamicVCC.Components.Units.HX.Piping pipe(
    redeclare final package Medium_1 = Medium_1,
    redeclare final model HeatTransfer_1 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.ConstHTC (alpha0=100),
    redeclare final model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.HeatLossAmbient (alpha0=10),
    redeclare final model Friction_1 = Friction_1_OD,
    final modelStructure=DynamicVCC.Components.Types.ModelStructure.av_vb,
    final EnableReverseFlow=false,
    final SteadyState_init=SteadyState_init,
    final useLumpedPressure=useLumpedPressure,
    final use_I_flows=use_I_flows,
    final Ncell=N_pipe,
    final T_amb=T_amb,
    final p_scale=p_scale,
    final h_scale=h_scale,
    final T_scale=T_scale,
    d_i=0.0110744,
    d_o=0.0127,
    L_tube=1.0922,
    p_init=fill(p_init_ID[1], N_pipe),
    h_init=fill(h_init_ID[1], N_pipe),
    Tt_init=fill(Tt_init_ID[1], N_pipe),
    m_flows_init=fill(m_flows_init_1[1], N_pipe + 1)) annotation (Placement(transformation(extent={{130,60},{110,80}})));

  DynamicVCC.Components.Units.MassFlowDevices.Accumulator accumulator(
    redeclare package Medium = Medium_1,
    m_flow_init=m_flows_init_1[1],
    V=0.0058,
    V_f_init=0.5*0.0058,
    p_init=p_init_OD[Ncell],
    SteadyState_init=SteadyState_init) annotation (Placement(transformation(extent={{108,-78},{142,-44}})));

  SI.Temperature T_amb;

  Modelica.Fluid.Sources.MassFlowSource_T fan_ID[Ncell](
    redeclare each final package Medium = Medium_2,
    each use_m_flow_in=true,
    each use_T_in=true,
    each use_X_in=true,
    each nPorts=1) annotation (Placement(transformation(extent={{-88,22},{-68,42}})));

  Modelica.Fluid.Sources.MassFlowSource_T fan_OD[Ncell](
    redeclare each final package Medium = Medium_2,
    each use_m_flow_in=true,
    each use_T_in=true,
    each use_X_in=true,
    each nPorts=1) annotation (Placement(transformation(extent={{-88,-110},{-68,-90}})));

  Modelica.Fluid.Sources.Boundary_pT sink_ID(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{-106,88},{-86,108}})));
  Modelica.Fluid.Sources.Boundary_pT sink_OD(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{-92,-24},{-72,-4}})));

  /* Sensor */
  DynamicVCC.Components.Units.Sensors.T_Superheat superheat(redeclare package Medium = Medium_1) annotation (Placement(transformation(extent={{40,-86},{60,-66}})));
  DynamicVCC.Components.Units.Sensors.T_Superheat subcooling(redeclare package Medium = Medium_1) annotation (Placement(transformation(extent={{-148,78},{-128,98}})));
  DynamicVCC.Components.Units.Sensors.T_Superheat superheat_suc(redeclare package Medium = Medium_1) "Suction superheat" annotation (Placement(transformation(extent={{176,-72},{196,-52}})));

  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply_ID(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{20,84},{40,104}})));
  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply_OD(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{-120,-42},{-100,-22}})));

/* Actuations */
  Modelica.Blocks.Sources.Constant ReversingValveOpen(k=1.0) annotation (Placement(transformation(extent={{136,-30},{116,-10}})));
  Modelica.Blocks.Continuous.FirstOrder CompSpeedFilter(
    k=1,
    T=1,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=speed_init) "First order filter of compressor speed" annotation (Placement(transformation(extent={{200,78},{180,98}})));
  Modelica.Blocks.Sources.Constant speed(k=53) "Input compressor speed" annotation (Placement(transformation(extent={{230,114},{210,134}})));

  Modelica.Blocks.Continuous.FirstOrder EXVOpeningFilter(
    T=1,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=opening_init) annotation (Placement(transformation(extent={{-206,-30},{-186,-10}})));
  Modelica.Blocks.Nonlinear.Limiter EXVOpeningLimiter(uMax=1.0, uMin=0.1) annotation (Placement(transformation(extent={{-184,12},{-164,32}})));
  Modelica.Blocks.Sources.Constant opening(k=0.16) "Input compressor speed" annotation (Placement(transformation(extent={{-236,12},{-216,32}})));

equation

  /* Control signals */
    //Comp_speed.u=if T_amb<280.4 then -2.086*T_amb+639.2 else 53.3;

    fan_ID_flow.u=0.00667*T_amb-1.5;

   for i in 1:Ncell loop
     fan_OD[i].m_flow_in=1.6/Ncell;
     fan_OD[i].X_in={0.0022,1-0.0022};
     fan_OD[i].T_in=293;
   end for;

     for i in 1:Ncell loop
   fan_ID[i].m_flow_in=fan_ID_flow.y/Ncell;
   fan_ID[i].T_in=295;
   fan_ID[i].X_in={0.002,1-0.002};
    end for;

   T_amb=293;

  /*
   connect(sampler_sh.u,superheat.T);
   exv_controller.u=1e-3-sampler_sh.y;
   connect(exv_controller.y,ZOH_exv.u);
   connect(ZOH_exv.y,exv_limiter.u);
   connect(exv_limiter.y,exv.opening);
*/

  connect(IndoorCoil.port_b1, exv.port_a) annotation (Line(points={{-56,66},{-134,66},{-134,36}}, color={0,127,255}));
  connect(exv.port_b, OutdoorCoil.port_a1) annotation (Line(points={{-134,8},{-134,-56},{-52,-56}},   color={0,127,255}));
  connect(fan_ID.ports[1], IndoorCoil.ports_a2) annotation (Line(points={{-68,32},{-68,12},{-20,12},{-20,41.52}}, color={0,127,255}));
  connect(sink_ID.ports, IndoorCoil.ports_b2) annotation (Line(points={{-86,98},{-62,98},{-62,108},{-20,108},{-20,90.48}}, color={0,127,255}));
  connect(sink_OD.ports, OutdoorCoil.ports_b2) annotation (Line(points={{-72,-14},{-16,-14},{-16,-31.52}}, color={0,127,255}));
  connect(superheat.port, OutdoorCoil.port_b1) annotation (Line(points={{50,-86},{50,-92},{28,-92},{28,-56},{20,-56}}, color={0,127,255}));
  connect(subcooling.port, IndoorCoil.port_b1) annotation (Line(points={{-138,78},{-138,66},{-56,66}}, color={0,127,255}));
  connect(T_supply_ID.ports, IndoorCoil.ports_b2) annotation (Line(points={{30,84.4},{30,76},{50,76},{50,110},{-20,110},{-20,90.48}}, color={0,127,255}));
  connect(T_supply_OD.ports, OutdoorCoil.ports_b2) annotation (Line(points={{-110,-41.6},{-110,-44},{-62,-44},{-62,-14},{-16,-14},{-16,-31.52}}, color={0,127,255}));
  connect(compressor.port_b, pipe.port_a1) annotation (Line(points={{150.364,39.6667},{150.364,70},{130,70}},
                                                                                                 color={0,127,255}));
  connect(pipe.port_b1, reversingValve.port_a) annotation (Line(points={{110,70},{102,70},{102,12},{79,12}},                    color={0,127,255}));
  connect(reversingValve.port_b2, IndoorCoil.port_a1) annotation (Line(points={{66.46,-26},{66.46,-32},{26,-32},{26,66},{16,66}}, color={0,127,255}));
  connect(OutdoorCoil.port_b1, reversingValve.port_a2) annotation (Line(points={{20,-56},{92,-56},{92,-26},{91.54,-26}}, color={0,127,255}));
  connect(reversingValve.port_b, accumulator.port_a) annotation (Line(points={{79,-26},{106,-26},{106,-38},{98,-38},{98,-61},{108,-61}}, color={0,127,255}));
  connect(accumulator.port_b, compressor.port_a) annotation (Line(points={{142,-61},{150.364,-61},{150.364,-3.66667}}, color={0,127,255}));
  connect(superheat_suc.port, accumulator.port_b) annotation (Line(points={{186,-72},{186,-78},{152,-78},{152,-61},{142,-61}}, color={0,127,255}));
  connect(fan_OD.ports[1], OutdoorCoil.ports_a2) annotation (Line(points={{-68,-100},{-16,-100},{-16,-80.48}}, color={0,127,255}));
  connect(ReversingValveOpen.y, reversingValve.opening_SP) annotation (Line(points={{115,-20},{106,-20},{106,2.88},{96.1,2.88}}, color={0,0,127}));
  connect(CompSpeedFilter.y, compressor.speed) annotation (Line(points={{179,88},{170,88},{170,50},{180,50},{180,36.85},{169.509,36.85}}, color={0,0,127}));
  connect(speed.y, CompSpeedFilter.u) annotation (Line(points={{209,124},{202,124},{202,104},{210,104},{210,88},{202,88}}, color={0,0,127}));
  connect(EXVOpeningLimiter.y, exv.opening) annotation (Line(points={{-163,22},{-140.44,21.86}}, color={0,0,127}));
  connect(EXVOpeningFilter.y, EXVOpeningLimiter.u) annotation (Line(points={{-185,-20},{-178,-20},{-178,6},{-194,6},{-194,22},{-186,22}}, color={0,0,127}));
  connect(opening.y, EXVOpeningFilter.u) annotation (Line(points={{-215,22},{-208,22},{-208,0},{-218,0},{-218,-20},{-208,-20}}, color={0,0,127}));
  annotation (experiment(
      StartTime=6000,
      StopTime=15000,
      Tolerance=0.001,
      __Dymola_Algorithm="Dassl"),
    Diagram(coordinateSystem(extent={{-240,-140},{240,160}})),
    Icon(coordinateSystem(extent={{-240,-140},{240,160}})));
end AirSourceHeatPump;