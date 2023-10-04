within DynamicVCC.Test.ModelCalibration;
model TestCycle "Cycle model for calibration"
  extends Modelica.Icons.Example;

  Modelica.Blocks.Interfaces.RealOutput y[7];

  // Components
  extends DynamicVCC.Test.ModelCalibration.HeatExchangers(
    final p_a_start_OD=p_a_start_OD_ss,
    final p_b_start_OD=p_b_start_OD_ss,
    final h_init_OD=h_init_OD_ss,
    final Tt_init_OD=Tt_init_OD_ss,
    final p_a_start_ID=p_a_start_ID_ss,
    final p_b_start_ID=p_b_start_ID_ss,
    final h_init_ID=h_init_ID_ss,
    final Tt_init_ID=Tt_init_ID_ss,
    final m_flow_init=m_flow_init_ss);

  // Initial conditions
  parameter Real speed_init=58.3;
  parameter Real opening_init=0.4821;

  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=2.3879e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=2.2985e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={445203.250000000,431486.250000000,418800.250000000,406256,393853.781250000,381593.781250000,369476.218750000,357501.343750000,345669.375000000,333980.562500000,322435.093750000,311033.187500000,299775.125000000,288661.093750000,277691.343750000,266866.062500000,260609.890625000,256855.140625000,254619.546875000,253294.906250000};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={319.076660156250,312.685058593750,312.107421875000,312.028015136719,311.948455810547,311.868774414063,311.788970947266,311.709014892578,311.628967285156,311.548767089844,311.468444824219,311.388000488281,311.307434082031,311.226715087891,311.145904541016,311.064941406250,308.505065917969,307.103637695313,306.252502441406,305.742156982422};

  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=1.0151e6;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=9.3448e5;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={262631.656250000,272160.500000000,281787.500000000,291512.968750000,301337.343750000,311260.968750000,321284.250000000,331407.531250000,341631.250000000,351955.781250000,362381.531250000,372908.906250000,383538.281250000,394270.062500000,405104.718750000,416042.625000000,422953.406250000,427742,431257.968750000,433838.843750000};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={281.732330322266,281.600006103516,281.467163085938,281.333801269531,281.199951171875,281.065582275391,280.930694580078,280.795257568359,280.659332275391,280.522827148438,280.385833740234,280.248260498047,280.110137939453,279.971466064453,279.832244873047,279.692474365234,285.144317626953,288.017303466797,289.740142822266,291.006042480469};

  parameter Medium_1.MassFlowRate m_flow_init_ss=0.052;

  DynamicVCC.Test.JCI.Compressor compressor(
  redeclare package Medium=Medium_1,
  Vs=2.9438902e-05,
  m_flow_init=system.m_flow_init,
  m_flow_nominal=0.05,
  speed_nominal=58.3,
  p_dis_init=p_a_start_OD,
  p_suc_init=p_b_start_ID,
  h_dis_init=h_init_OD[1],
  h_suc_init=h_init_ID[Ncell]) annotation (Placement(transformation(extent={{112,28},{70,70}})));

  DynamicVCC.Test.JCI.EXV exv(
  redeclare package Medium=Medium_1,
  Av=3.1416e-6,
  p_a_init=p_b_start_OD,
  p_b_init=p_a_start_ID,
  final dp_nominal=10e5) annotation (Placement(transformation(extent={{14,14},{-14,-14}},
        rotation=90,
        origin={-80,-40})));

  Modelica.Fluid.Sources.MassFlowSource_T airsource_OD[Ncell](
  redeclare each final package Medium=Medium_2,
  each use_m_flow_in=true,
  each use_T_in=true,
  each use_X_in=true,
  each nPorts=1) annotation (Placement(transformation(extent={{12,0},{-12,24}})));

  Modelica.Fluid.Sources.Boundary_pT airsink_OD[Ncell](
  redeclare each package Medium=Medium_2,
  each nPorts=1) annotation (Placement(transformation(extent={{44,70},{26,88}})));

  Modelica.Fluid.Sources.MassFlowSource_T airsource_ID[Ncell](
  redeclare each final package Medium=Medium_2,
  each use_m_flow_in=true,
  each use_T_in=true,
  each use_X_in=true,
  each nPorts=1) annotation (Placement(transformation(extent={{72,-98},{52,-74}})));

  Modelica.Fluid.Sources.Boundary_pT airsink_ID[Ncell](
  redeclare each package Medium=Medium_2,
  each nPorts=1) annotation (Placement(transformation(extent={{20,-32},{2,-14}})));

  DynamicVCC.Components.Units.Sensors.T_Superheat superheat(
  redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{144,24},{122,46}})));

  DynamicVCC.Components.Units.Sensors.T_Superheat subcooling(
  redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-72,66},{-94,88}})));

  // Connect piping between indoor and outdoor units
  replaceable model HeatTransfer_1_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    alpha0=50);

  replaceable model HeatTransfer_2_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
   final alpha0=5) "Connect piping heat loss to ambient";

  replaceable model FlowModel_1_LL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
  final dp_nominal=1.3796e4/Ncell_piping,
  final k=1,
  final b=2);

  replaceable model FlowModel_1_VL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
  final dp_nominal=2.0953e4/Ncell_piping,
  final k=1,
  final b=2);

  parameter Integer Ncell_piping=3;

  Components.Units.HX.Piping liquidLine(
    redeclare final package Medium_1=Medium_1,
    final Ncell=Ncell_piping,
    redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
    redeclare final model FlowModel_1=FlowModel_1_LL,
    redeclare final model SlipRatio=SlipRatio,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final useLumpedPressure=useLumpedPressure,
    d_o=0.009525,
    d_i=0.0078994,
    length=6.1,
    T_amb=300,
    p_a_start=p_b_start_OD,
    h_init=fill(h_init_OD[Ncell],Ncell_piping),
    Tt_init=fill(Tt_init_OD[Ncell], Ncell_piping))
                         annotation (Placement(transformation(
        extent={{12,-12},{-12,12}},
        rotation=90,
        origin={-80,38})));

  Components.Units.HX.Piping vaporLine(
  redeclare final package Medium_1=Medium_1,
    final Ncell=Ncell_piping,
    redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
    redeclare final model FlowModel_1=FlowModel_1_VL,
    redeclare final model SlipRatio=SlipRatio,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final useLumpedPressure=useLumpedPressure,
    d_o=0.009525,
    d_i=0.0078994,
    length=6.1,
    T_amb=300,
    p_a_start=p_b_start_ID,
    h_init=fill(h_init_ID[Ncell],Ncell_piping),
    Tt_init=fill(Tt_init_ID[Ncell], Ncell_piping))       annotation (Placement(transformation(
        extent={{12,-12},{-12,12}},
        rotation=270,
        origin={112,-26})));

/*
  Components.Units.MassFlowDevices.Valve.ReversingValve reversingValve(
    redeclare final package Medium=Medium_1,
    hpMode=DynamicVCC.Components.Types.HPmode.heating,
    p_dis_init=p_a_start_OD,
    p_suc_init=p_b_start_ID,
    A_open=1.0e-4,
    dp_nominal=20000,
    C_Hd=1,
    C_Hs=0.7) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={52,22})));

  Components.Units.HX.Piping suctionLine(
    redeclare final package Medium_1 = Medium_1,
    final Ncell=Ncell_piping,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_piping,
    redeclare final model FlowModel_1 = FlowModel_1_VL,
    redeclare final model SlipRatio = SlipRatio,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final useLumpedPressure=useLumpedPressure,
    d_o=0.009525,
    d_i=0.0078994,
    length=0.1,
    T_amb=300,
    p_a_start=p_b_start_ID,
    h_init=fill(h_init_ID[Ncell], Ncell_piping),
    Tt_init=fill(Tt_init_ID[Ncell], Ncell_piping)) annotation (Placement(transformation(
        extent={{8,-8},{-8,8}},
        rotation=180,
        origin={70,-6})));

  Components.Units.HX.Piping dischargeLine(
    redeclare final package Medium_1 = Medium_1,
    final Ncell=Ncell_piping,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_piping,
    redeclare final model FlowModel_1 = FlowModel_1_VL,
    redeclare final model SlipRatio = SlipRatio,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final useLumpedPressure=useLumpedPressure,
    d_o=0.009525,
    d_i=0.0078994,
    length=0.1,
    T_amb=300,
    p_a_start=p_a_start_OD,
    h_init=fill(h_init_OD[1], Ncell_piping),
    Tt_init=fill(Tt_init_OD[1], Ncell_piping)) annotation (Placement(transformation(
        extent={{-7,-7},{7,7}},
        rotation=180,
        origin={49,49})));
     

  Modelica.Fluid.Sources.Boundary_ph boundary1(
    redeclare package Medium = Medium_1,
    use_p_in=true,
    nPorts=1) annotation (Placement(transformation(extent={{-116,6},{-100,22}})));

  Modelica.Fluid.Sources.Boundary_ph boundary2(
    redeclare package Medium = Medium_1,
    use_p_in=true,
    use_h_in=true,
    nPorts=1) annotation (Placement(transformation(extent={{-116,-26},{-100,-10}})));
 */
  Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(extent={{-128,-12},{-108,8}})));
  Modelica.Blocks.Sources.Constant const(k=10) annotation (Placement(transformation(extent={{-164,-12},{-144,8}})));
  Modelica.Blocks.Continuous.PI PI(
    k=-0.001,
    T=5,
    initType=Modelica.Blocks.Types.Init.InitialState,
    x_start=-482.129119873047)
                          annotation (Placement(transformation(extent={{-144,-76},{-124,-56}})));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax=1.0, uMin=0.01) annotation (Placement(transformation(extent={{-116,-76},{-96,-56}})));
  Modelica.Blocks.Sources.Pulse pulse(
    amplitude=-15,
    width=50,
    period=1000,
    offset=58.3,
    startTime=300) annotation (Placement(transformation(extent={{38,22},{58,42}})));

  SI.Mass charge;

   Modelica.Blocks.Sources.CombiTimeTable OutputData(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,tableName="outputData",fileName="C:/Jiacheng Ma/outputData.mat",columns=2:8);

equation

  charge=IndoorCoil.charge+OutdoorCoil.charge+liquidLine.charge+vaporLine.charge;
  compressor.T_amb=airsource_OD[1].T_in;

  y[1]=compressor.port_b.p;
  y[2]=compressor.port_a.p;
  y[3]=IndoorCoil.Ta_out_ave;
  y[4]=OutdoorCoil.Ta_out_ave;
  y[5]=subcooling.T;
  y[6]=compressor.Pwr;
  y[7]=IndoorCoil.Q_flow_2;

  //openingFilter.u=BC.y[2];
  //compressorSpeedFilter.u=BC.y[1];

  //firstOrder.u=25;//BC.y[1];
  //compressor.speed=u[1];//BC.y[1];
  //exv.opening=u[2];//BC.y[2];

  for i in 1:Ncell loop
    airsource_ID[i].m_flow_in=1/Ncell;
    airsource_ID[i].T_in=294.5;
    airsource_ID[i].X_in={1e-5,1-1e-5};
    airsource_OD[i].m_flow_in=2.1/Ncell;
    airsource_OD[i].T_in=305;
    airsource_OD[i].X_in={1e-5,1-1e-5};
  end for;

  connect(OutdoorCoil.ports_a2,airsource_OD.ports[1]);
  connect(OutdoorCoil.ports_b2,airsink_OD.ports[1]);
  connect(IndoorCoil.ports_a2,airsource_ID.ports[1]);
  connect(IndoorCoil.ports_b2,airsink_ID.ports[1]);
  connect(subcooling.port, OutdoorCoil.port_b1) annotation (Line(points={{-83,66},{-83,60},{-22,60}},                   color={0,127,255},
      thickness=0.5));
  connect(compressor.port_b, OutdoorCoil.port_a1) annotation (Line(points={{70,49},{70,48},{30,48},{30,60},{22,60}}, color={0,127,255}));
  connect(exv.port_b, IndoorCoil.port_a1) annotation (Line(points={{-80,-54},{-80,-60},{-22,-60}}, color={0,127,255}));
  connect(superheat.port, compressor.port_a) annotation (Line(points={{133,24},{134,24},{134,16},{112,16},{112,22},{118,22},{118,49},{112,49}}, color={0,127,255}));
  connect(exv.port_a, liquidLine.port_b1) annotation (Line(points={{-80,-26},{-80,26}}, color={0,127,255}));
  connect(liquidLine.port_a1, OutdoorCoil.port_b1) annotation (Line(points={{-80,50},{-80,60},{-22,60}}, color={0,127,255}));
  connect(IndoorCoil.port_b1, vaporLine.port_a1) annotation (Line(points={{22,-60},{112,-60},{112,-38}}, color={0,127,255}));
  connect(vaporLine.port_b1, compressor.port_a) annotation (Line(points={{112,-14},{114,-14},{114,16},{112,16},{112,22},{118,22},{118,49},{112,49}}, color={0,127,255}));
  connect(superheat.T, feedback.u2) annotation (Line(points={{125.3,35},{124,35},{124,12},{22,12},{22,-8},{-102,-8},{-102,-18},{-118,-18},{-118,-10}}, color={0,0,127}));
  connect(const.y, feedback.u1) annotation (Line(points={{-143,-2},{-126,-2}}, color={0,0,127}));
  connect(feedback.y, PI.u) annotation (Line(points={{-109,-2},{-106,-2},{-106,-44},{-150,-44},{-150,-66},{-146,-66}}, color={0,0,127}));
  connect(PI.y, limiter.u) annotation (Line(points={{-123,-66},{-118,-66}}, color={0,0,127}));
  connect(limiter.y, exv.opening) annotation (Line(points={{-95,-66},{-60,-66},{-60,-20},{-100,-20},{-100,-40.14},{-86.44,-40.14}}, color={0,0,127}));
  connect(pulse.y, compressor.speed) annotation (Line(points={{59,32},{59,34.825},{72.73,34.825}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(extent={{-180,-100},{160,100}})), Icon(coordinateSystem(extent={{-180,-100},{160,100}})),
    experiment(
      StartTime=1,
      StopTime=1200,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end TestCycle;
