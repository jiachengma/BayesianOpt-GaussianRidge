within DynamicVCC.Test.JCI.Test.TPWL;
model TestCycle_LoadChange
  extends Modelica.Icons.Example;

  input Real u[6](start={25,0.23,1,2.1,294.5,305});

  //Modelica.Blocks.Interfaces.RealOutput charge annotation (Placement(transformation(extent={{160,70},{132,92}}), iconTransformation(extent={{160,70},{132,92}})));
  Modelica.Blocks.Interfaces.RealOutput y[8];
  //Modelica.Blocks.Interfaces.RealOutput COP annotation (Placement(transformation(extent={{160,50},{132,72}}), iconTransformation(extent={{160,50},{132,72}})));

  // Components
  extends DynamicVCC.Test.JCI.Test.TPWL.HeatExchangers(
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
  parameter Real opening_init=0.458;

  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=2.4144e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=2.3312e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={444721.687500000,431018.437500000,418197.875000000,405618.281250000,393263.406250000,381117.125000000,369163.218750000,357385.625000000,345768.187500000,334294.875000000,322949.625000000,311716.406250000,300579.250000000,289522.125000000,278529.062500000,267584.062500000,261193.046875000,257376.937500000,255116.484375000,253783.609375000};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={319.365295410156,312.924743652344,312.433593750000,312.299530029297,312.174468994141,312.058410644531,311.951354980469,311.853271484375,311.764160156250,311.683959960938,311.612701416016,311.550354003906,311.496917724609,311.452392578125,311.416748046875,311.389984130859,308.856079101563,307.423339843750,306.557769775391,306.041625976563};

  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=9.6703e5;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=9.4187e5;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={266489.250000000,277491.781250000,288509.781250000,299545.500000000,310601.125000000,321678.937500000,332781.218750000,343910.187500000,355068.218750000,366257.562500000,377480.625000000,388739.718750000,400037.281250000,411375.718750000,422757.468750000,427689.718750000,431256.906250000,433832.937500000,435692.843750000,437037.375000000};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={280.202697753906,280.185455322266,280.165313720703,280.142272949219,280.116333007813,280.087493896484,280.055694580078,280.020904541016,279.983154296875,279.942382812500,279.898529052734,279.851623535156,279.801605224609,279.748413085938,279.692047119141,288.083007812500,289.859039306641,291.148529052734,292.080200195313,292.750701904297};
  parameter Medium_1.MassFlowRate m_flow_init_ss=0.05;

/*
  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=2.1647e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=2.1344e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={438136.125000000,427509.187500000,414708.687500000,402018.968750000,389440.125000000,376972.250000000,364615.375000000,352369.468750000,340234.593750000,328210.687500000,316297.812500000,304495.875000000,292804.781250000,281224.281250000,269753.593750000,258390.218750000,254478.750000000,252796.265625000,252076.234375000,251769.265625000};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={311.998138427734,307.689361572266,308.238098144531,308.209991455078,308.181823730469,308.153686523438,308.125518798828,308.097320556641,308.069122314453,308.040924072266,308.012695312500,307.984436035156,307.956176757813,307.927917480469,307.899627685547,307.871368408203,305.986785888672,305.424224853516,305.181457519531,305.077362060547};

  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=1.0527e6;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=1.0391e6;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={270200.062500000,288846.437500000,307530.906250000,326252.843750000,345011.687500000,363807.187500000,382639.031250000,401507.031250000,420411.062500000,427479.906250000,431733.937500000,434283.406250000,435809.125000000,436723.125000000,437272.687500000,437605.562500000,437809.750000000,437937.562500000,438020.093750000,438075.750000000};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={283.094940185547,283.073730468750,283.052520751953,283.031311035156,283.010101318359,282.988830566406,282.967590332031,282.946319580078,282.925048828125,290.250823974609,291.943572998047,292.968322753906,293.583648681641,293.951232910156,294.170166015625,294.300292968750,294.377563476563,294.423400878906,294.450592041016,294.466735839844};
  parameter Medium_1.MassFlowRate m_flow_init_ss=0.0235;
*/

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

  Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    tableName="BC",
    fileName="C:/Jiacheng Ma/BoundaryCondition/JCI/4_30_2021/BC.mat",                                                                                                                                     columns=2:6) annotation (Placement(transformation(extent={{-124,36},{-146,58}})));
  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,                                                                tableName="Mea",
    fileName="C:/Jiacheng Ma/BoundaryCondition/JCI/4_30_2021/Mea.mat",                                                                                                                                 columns=2:10) annotation (Placement(transformation(extent={{-124,70},{-146,92}})));

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
  /*
  replaceable model FlowModel_1_LL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  lambda0=0.1);

  replaceable model FlowModel_1_VL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  lambda0=0.1);
 */
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


    SI.MassFlowRate m_flow_air_ID;
    SI.MassFlowRate m_flow_air_OD;
    SI.Temperature T_ID;
    SI.Temperature T_OD;
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
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=opening_init) annotation (Placement(transformation(extent={{-144,-76},{-124,-56}})));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax=1.0, uMin=0.01) annotation (Placement(transformation(extent={{-116,-76},{-96,-56}})));
  /*
  Modelica.Blocks.Sources.Pulse pulse(
    amplitude=-15,
    width=50,
    period=1000,
    offset=58.3,
    startTime=5000)
                   annotation (Placement(transformation(extent={{38,22},{58,42}})));
*/
  Real charge;
  Real COP;

initial equation
  //charge=4;

equation

  charge=IndoorCoil.charge+OutdoorCoil.charge+liquidLine.charge+vaporLine.charge;
  COP=abs(IndoorCoil.Q_flow_2)/(compressor.Pwr+1e-5);
  compressor.T_amb=airsource_OD[1].T_in;

  y[1]=compressor.port_b.p;
  y[2]=compressor.port_a.p;
  y[3]=IndoorCoil.Ta_out_ave;
  y[4]=OutdoorCoil.Ta_out_ave;
  y[5]=superheat.T;
  y[6]=subcooling.T;
  y[7]=compressor.Pwr;
  y[8]=IndoorCoil.Q_flow_2;

  //openingFilter.u=BC.y[2];
  //compressorSpeedFilter.u=BC.y[1];
  compressor.speed=58.3;

  //firstOrder.u=25;//BC.y[1];
  //compressor.speed=u[1];//BC.y[1];
  //exv.opening=u[2];//BC.y[2];

/*
  boundary2.p_in=21.0018e5;
  boundary2.h_in=251595;
  boundary1.p_in=21.0018e5;
*/

  for i in 1:Ncell loop
    airsource_ID[i].m_flow_in=m_flow_air_ID/Ncell;
    airsource_ID[i].T_in=T_ID;
    airsource_ID[i].X_in={1e-4,1-1e-4};
    airsource_OD[i].m_flow_in=m_flow_air_OD/Ncell;
    airsource_OD[i].T_in=T_OD;
    airsource_OD[i].X_in={1e-4,1-1e-4};
  end for;

  m_flow_air_ID=u[3];
  m_flow_air_OD=u[4];
  T_ID=u[5];
  T_OD=u[6];

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
  annotation (Diagram(coordinateSystem(extent={{-180,-100},{160,100}})), Icon(coordinateSystem(extent={{-180,-100},{160,100}})),
    experiment(
      StartTime=1,
      StopTime=2000,
      __Dymola_Algorithm="Dassl"));
end TestCycle_LoadChange;
