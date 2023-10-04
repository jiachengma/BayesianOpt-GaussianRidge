within DynamicVCC.Examples.Tests;
model CompressorOnOffCycle "Test compressor on/off cycles"
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

  parameter SI.CoefficientOfHeatTransfer alpha_nominal_ID=8000;
  parameter SI.CoefficientOfHeatTransfer alpha_nominal_OD=5000;
  parameter SI.LewisNumber Le=0.854 "Lewis number";
  parameter Integer N_pipe=5; //Discharge pipe line
  parameter SI.Time samplePeriod_exv=90;
  parameter SI.Time startTime=0 "Simulation start time";

/***************************** Initialization *******************************/
/*
      parameter Real p_init_ID[Ncell]=linspace(33.6e5,33.5e5,Ncell);
      parameter Real h_init_ID[Ncell]=linspace(4.35e5,2.5e5,Ncell);
      parameter Real Tt_init_ID[Ncell]=linspace(318,302,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_ID=0.435;
      parameter SI.Temperature T_a_init_ID[Ncell]=linspace(317,301,Ncell);
      parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.002,Ncell);

      parameter Real p_init_OD[Ncell]=linspace(12.7e5,11.7e5,Ncell);
      parameter Real h_init_OD[Ncell]=linspace(2.58e5,4.25e5,Ncell);
      parameter Real Tt_init_OD[Ncell]=linspace(288.7,287.7,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_OD=1.6;
      parameter SI.Temperature T_a_init_OD[Ncell]=linspace(288,287,Ncell);
      parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.002,Ncell);
      */

      //parameter Real p_init_ID[Ncell]=fill(3.2e6,Ncell);
      parameter Real p_init_ID[Ncell]={3259083.98628235,3257914.06631470,3257197.14164734,3256520.74813843,3255886.79313660,3255295.75347900,3254746.91390991,3254240.98968506,3253777.50396729,3253356.45675659,3252978.08647156,3252642.15469360,3252348.89984131,3252097.84507751,3251889.22882080,3251722.81265259,3251598.59657288,3251475.57258606,3251354.93278503,3251236.43875122,3251119.85206604,3251004.69589233,3250891.20864868,3250778.67507935,3250667.33360291,3250556.70738220,3250447.03483582,3250338.07754517,3250229.59709168,3250067.71087647};
      parameter Real h_init_ID[Ncell]={434678.888320923,428404.569625855,419242.668151855,410085.868835449,400933.980941772,391786.813735962,382644.200325012,373505.973815918,364371.967315674,355241.990089417,346115.875244141,336993.479728699,327874.588966370,318759.059906006,309646.701812744,300537.347793579,291430.830955505,285892.415046692,280869.913101196,276316.165924072,272189.569473267,268452.405929565,265070.128440857,262011.051177979,259245.920181274,256747.889518738,254492.330551147,252456.641197205,250620.150566101,248964.047431946};
      parameter Real Tt_init_ID[Ncell]={316.273612260819,315.033094882965,324.252750992775,324.236511468887,324.220857739449,324.205757260323,324.191210031509,324.177183508873,324.163742780685,324.150822758675,324.138488531113,324.126675009728,324.115414738655,324.104707717896,324.094553947449,324.084953427315,324.075873613358,312.683701157570,311.036414623261,309.539709806442,308.175882697105,306.932468891144,305.799249529839,304.767405152321,303.828832268715,302.975980639458,302.201820731163,301.499778628349,300.863670945168,300.287802457809};
      parameter SI.MassFlowRate m_dot_air_ini_ID=0.435;
      parameter SI.Temperature T_a_init_ID[Ncell]={315.658325195313,314.453704833984,323.406677246094,323.390930175781,323.375732421875,323.361053466797,323.346923828125,323.333312988281,323.320251464844,323.307739257813,323.295745849609,323.284271240234,323.273345947266,323.262939453125,323.253082275391,323.243743896484,323.234954833984,312.172271728516,310.572601318359,309.119201660156,307.794830322266,306.587341308594,305.486907958984,304.484893798828,303.573455810547,302.745269775391,301.993530273438,301.311767578125,300.694091796875,300.134857177734};
      parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.002,Ncell);

      //parameter Real p_init_OD[Ncell]=fill(12e5,Ncell);
      parameter Real p_init_OD[Ncell]={1258797.04952240,1257513.52310181,1256555.55725098,1255498.88610840,1254342.43679047,1253085.37483215,1251726.50814056,1250264.64462280,1248698.71139526,1247027.15873718,1245248.67534637,1243361.83071136,1241364.83669281,1239256.02436066,1237033.72478485,1234695.91140747,1232240.55767059,1229665.51780701,1226968.52684021,1224147.20058441,1221199.03564453,1218121.40941620,1214911.34166718,1211565.97137451,1208082.19909668,1204456.68697357,1200686.09714508,1196766.73412323,1192694.66400146,1186402.08244324};
      parameter Real h_init_OD[Ncell]={254098.701477051,259268.379211426,264463.448524475,269686.818122864,274941.468238831,280230.379104614,285556.626319885,290923.237800598,296333.360671997,301790.189743042,307296.872138977,312856.745719910,318473.124504089,324149.370193481,329888.963699341,335695.457458496,341572.403907776,347523.498535156,353552.532196045,359663.343429565,365859.889984131,372146.201133728,378526.449203491,385004.901885986,391585.922241211,398274.087905884,405073.976516724,411990.451812744,419028.377532959,424708.318710327};
      parameter Real Tt_init_OD[Ncell]={288.456153988838,288.425171971321,288.402716517448,288.377690076828,288.350027561188,288.319696426392,288.286696672440,288.250963211060,288.212463498116,288.171164989471,288.127035140991,288.080008864403,288.030021071434,287.977006673813,287.920965671539,287.861800432205,287.799445867538,287.733836889267,287.664875864983,287.592497706413,287.516669869423,287.437229633331,287.354111909866,287.267219066620,287.176420927048,287.081652402878,286.982750773430,286.879618406296,286.772125124931,287.978601336479};
      parameter SI.MassFlowRate m_dot_air_ini_OD=1.6;
      parameter SI.Temperature T_a_init_OD[Ncell]={288.479858398438,288.449066162109,288.426727294922,288.401824951172,288.374298095703,288.344146728516,288.311309814453,288.275756835938,288.237487792969,288.196380615234,288.152465820313,288.105682373047,288.055969238281,288.003234863281,287.947479248047,287.888610839844,287.826599121094,287.761322021484,287.692718505859,287.620758056641,287.545288085938,287.466247558594,287.383575439453,287.297149658203,287.206817626953,287.112548828125,287.014190673828,286.911590576172,286.804626464844,288.004821777344};
      parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.002,Ncell);

      parameter SI.Frequency speed_init=53;
      parameter SI.Power Pwr_init=1900;
      parameter Real opening_init=0.16;
      parameter SI.Pressure dp_init=20e5;
      parameter SI.MassFlowRate m_flows_init_1[Ncell+1]=fill(0.047,Ncell+1);
      parameter SI.Volume V_f_init=0.0023;

/***************************** Numerical *******************************/
  import DynamicVCC.Components.Types.ModelStructure;
  import DynamicVCC.Components.Types.DifferentialState;

  parameter DifferentialState differentialState=DifferentialState.ph;
  parameter Boolean SteadyState_init=false;
  parameter Boolean EnableReverseFlow=false;
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
  alpha0=alpha_nominal_ID,
  alpha_cst=alpha_nominal_ID);

  replaceable model Evaporation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=alpha_nominal_OD,
  alpha_cst=alpha_nominal_ID);
 /*
  replaceable model SinglePhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.SinglePhase_Gnielinski (
  alpha0=alpha_nominal_OD);
*/
  replaceable model SinglePhase =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.Constant (
  alpha0=600,
  alpha_cst=600);

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
    final alpha0=5);

 replaceable model Friction_1_OD=DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant (       f0=0.1,dp_nominal=1e3,m_flow_nominal=0.048);

 replaceable model Friction_2_OD = DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
 f0=0.12,m_flow_nominal=1)
 "Outdoor air side pressure drop";

  replaceable model HeatTransfer_1_ID=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=SinglePhase,
  redeclare model VaporZone=SinglePhase,
  redeclare model TwoPhaseZone=Condensation);

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
    final alpha0=5);

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
 f0=0.1,m_flow_nominal=1)
 "Indoor air side pressure drop";

 //replaceable model VoidFraction = DynamicVCC.Components.Pipes.BaseClasses.VoidFraction.Zivi;

/********************************* Control ******************************************/

  Modelica.Blocks.Continuous.FirstOrder fan_ID_flow(T=15,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=m_dot_air_ini_ID);
  Modelica.Blocks.Continuous.FirstOrder fan_OD_flow(T=15,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=m_dot_air_ini_OD);

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
    m_flow_nominal=0.047,
    speed_nominal=53,
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
        origin={116,26})));

  DynamicVCC.Components.Units.MassFlowDevices.Valve.ElectronicExpansionValve exv(
    redeclare package Medium = Medium_1,
    final Av=2.5447e-06,
    final EnableReverseFlow=EnableReverseFlow,
    opening_init=opening_init,
    dp_nominal=20e5,
    m_flow_nominal=0.047,
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
    redeclare final model FreeConvection = FreeConvection_ID,
    redeclare final model Friction_1 = Friction_1_ID,
    redeclare final model Friction_2 = Friction_2_ID,
    final differentialState=differentialState,
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
    w_a_init=w_a_init_ID,
    SteadyState_init=SteadyState_init,
    m_flows_a_init=fill(m_dot_air_ini_ID/Ncell,Ncell),
    Ta_init=T_a_init_ID) annotation (Placement(transformation(extent={{16,30},{-56,102}})));

    DynamicVCC.Components.Units.HX.FinTubeHX OutdoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_OD,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_OD,
    redeclare final model FreeConvection = FreeConvection_OD,
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
    m_flows_a_init=fill(m_dot_air_ini_OD/Ncell,Ncell),
    Ta_init=T_a_init_OD) annotation (Placement(transformation(extent={{-52,-92},{20,-20}})));
/*
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

  DynamicVCC.Components.Pipes.ConnectingPipe pipe(
    redeclare final package Medium_1 = Medium_1,
    redeclare final model HeatTransfer_1 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstHTC (alpha0=100),
    redeclare final model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.HeatLossAmbient (alpha0=10),
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
    */

  DynamicVCC.Components.Units.MassFlowDevices.Accumulator accumulator(
    redeclare package Medium = Medium_1,
    m_flow_init=m_flows_init_1[1],
    V=0.0058,
    V_f_init=V_f_init,
    p_init=p_init_OD[Ncell],
    SteadyState_init=SteadyState_init) annotation (Placement(transformation(extent={{72,-52},{106,-18}})));

  SI.Temperature T_amb;

  SI.Mass charge;
  SI.Energy energy "Energy conserved in the cycle";

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

  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply_ID(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{20,84},{40,104}})));
  DynamicVCC.Components.Units.Sensors.Temperature_grid T_supply_OD(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{-120,-42},{-100,-22}})));

/* Actuations */
  //Modelica.Blocks.Sources.Constant ReversingValveOpen(k=1.0) annotation (Placement(transformation(extent={{136,-30},{116,-10}})));
  Modelica.Blocks.Continuous.FirstOrder CompSpeedFilter(
    k=1,
    T=20,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=speed_init) "First order filter of compressor speed" annotation (Placement(transformation(extent={{168,86},{148,106}})));

  Modelica.Blocks.Continuous.FirstOrder EXVOpeningFilter(
    T=1,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=opening_init) annotation (Placement(transformation(extent={{-206,-30},{-186,-10}})));
  Modelica.Blocks.Nonlinear.Limiter EXVOpeningLimiter(uMax=1.0, uMin=0.1) annotation (Placement(transformation(extent={{-184,12},{-164,32}})));
  Modelica.Blocks.Sources.Constant opening(k=0.16) "Input compressor speed" annotation (Placement(transformation(extent={{-236,12},{-216,32}})));

  Modelica.Blocks.Sources.Pulse pulse(
    amplitude=-53,
    width=50,
    period=2000,
    nperiod=5,
    offset=53.001,
    startTime=1000)
                   annotation (Placement(transformation(extent={{182,52},{202,72}})));
  Modelica.Blocks.Sources.Pulse pulse_fan(
    amplitude=-0.99,
    width=42,
    period=2000,
    nperiod=5,
    offset=1,
    startTime=1080) annotation (Placement(transformation(extent={{180,-60},{200,-40}})));
equation

  charge=IndoorCoil.charge+OutdoorCoil.charge+accumulator.charge;
  energy=IndoorCoil.energy+OutdoorCoil.energy+accumulator.energy;
  /* Control signals */

    CompSpeedFilter.u=pulse.y;//if T_amb<280 then -2.086*T_amb+639.2 else 53;
    //CompSpeedFilter.u=if time<1000 then 53 elseif time<2000 then 0.001 else 53;

    //fan_ID_flow.u=if time<1060 then 0.4543 elseif time<1940 then 1e-3 else 0.4543; //0.00667*T_amb-1.5;
    //fan_OD_flow.u=if time<1060 then 1.6 elseif time<1940 then 1e-3 else 1.6;

    fan_ID_flow.u=0.4543*pulse_fan.y;
    fan_OD_flow.u=1.6*pulse_fan.y;

   for i in 1:Ncell loop
     fan_OD[i].m_flow_in=fan_OD_flow.y/Ncell;
     fan_OD[i].X_in={0.0022,1-0.0022};
     fan_OD[i].T_in=T_amb;
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
  connect(fan_OD.ports[1], OutdoorCoil.ports_a2) annotation (Line(points={{-68,-100},{-16,-100},{-16,-80.48}}, color={0,127,255}));
  connect(CompSpeedFilter.y, compressor.speed) annotation (Line(points={{147,96},{138,96},{138,58},{148,58},{148,44.85},{137.509,44.85}}, color={0,0,127}));
  connect(EXVOpeningLimiter.y, exv.opening) annotation (Line(points={{-163,22},{-140.44,21.86}}, color={0,0,127}));
  connect(EXVOpeningFilter.y, EXVOpeningLimiter.u) annotation (Line(points={{-185,-20},{-178,-20},{-178,6},{-194,6},{-194,22},{-186,22}}, color={0,0,127}));
  connect(opening.y, EXVOpeningFilter.u) annotation (Line(points={{-215,22},{-208,22},{-208,0},{-218,0},{-218,-20},{-208,-20}}, color={0,0,127}));
  connect(compressor.port_b, IndoorCoil.port_a1) annotation (Line(points={{118.364,47.6667},{118.364,66},{16,66}}, color={0,127,255}));
  connect(OutdoorCoil.port_b1, accumulator.port_a) annotation (Line(points={{20,-56},{62,-56},{62,-35},{72,-35}}, color={0,127,255}));
  connect(accumulator.port_b, compressor.port_a) annotation (Line(points={{106,-35},{118.364,-35},{118.364,4.33333}}, color={0,127,255}));
  annotation (experiment(
      StopTime=12000,
      Interval=0.5,
      Tolerance=0.001,
      __Dymola_Algorithm="Dassl"),
    Diagram(coordinateSystem(extent={{-240,-140},{240,160}})),
    Icon(coordinateSystem(extent={{-240,-140},{240,160}})));
end CompressorOnOffCycle;
