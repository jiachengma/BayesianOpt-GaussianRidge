within DynamicVCC.Test.JCI.Test.TPWL;
model TestCycle_OnOff
  extends Modelica.Icons.Example;

  //input Real u[6](start={0,0.458,0,0,294.5,305});
  input Real u[6](start={58.3,0.458,1,2.1,294.5,305});

  Modelica.Blocks.Interfaces.RealOutput y[8];

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

  // steady-state
  parameter Real speed_init=58.3;
  parameter Real opening_init=0.458;
  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=2.4726e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=2.3841e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={447908.937500000,433267.531250000,419473,405802.125000000,392255.093750000,378832.125000000,365533.437500000,352359.218750000,339309.687500000,326385.031250000,313585.468750000,300911.218750000,288362.437500000,275939.343750000,265681.187500000,260356.390625000,256980.968750000,254855.312500000,253522.328125000,252688.906250000};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={320.383148193359,313.944549560547,313.427185058594,313.351654052734,313.275970458984,313.200195312500,313.124267578125,313.048217773438,312.972045898438,312.895751953125,312.819335937500,312.742797851563,312.666137695313,312.589355468750,311.266754150391,308.252960205078,307.062072753906,306.298583984375,305.814331054688,305.509124755859};

  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=1.0105e6;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=9.3031e5;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={262267.281250000,272034.843750000,281901.750000000,291868.375000000,301935.125000000,312102.375000000,322370.500000000,332739.875000000,343210.968750000,353784.093750000,364459.718750000,375238.187500000,386119.968750000,397105.468750000,408195.062500000,419389.218750000,425198.375000000,429450.156250000,432558.312500000,434831.406250000};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={281.588562011719,281.456390380859,281.323730468750,281.190551757813,281.056854248047,280.922668457031,280.787963867188,280.652709960938,280.516937255859,280.380615234375,280.243774414063,280.106414794922,279.968475341797,279.829986572266,279.690948486328,279.551330566406,286.742431640625,288.822174072266,290.349395751953,291.464508056641};

  parameter Medium_1.MassFlowRate m_flow_init_ss=0.0513;

  parameter Medium_1.AbsolutePressure p_a_start_LL=2.3841e6;
  parameter Medium_1.AbsolutePressure p_b_start_LL=2.3735e6;
  parameter Medium_1.SpecificEnthalpy h_init_LL[Ncell_piping]={252658.687500000,252628.625000000,252598.687500000};
  parameter SI.Temperature Tt_init_LL[Ncell_piping]={305.145202636719,305.128204345703,305.111297607422};

  parameter Medium_1.AbsolutePressure p_a_start_VL=9.3031e5;
  parameter Medium_1.AbsolutePressure p_b_start_VL=9.1419e5;
  parameter Medium_1.SpecificEnthalpy h_init_VL[Ncell_piping]={434889.937500000,434949.218750000,435009.25};
  parameter SI.Temperature Tt_init_VL[Ncell_piping]={290.209289550781,290.093780517578,289.978820800781};



 /*
  // off-cycle
  parameter Real speed_init=0;
  parameter Real opening_init=0.458;
  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=2.0534e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=2.0506e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={438158.215104701,431640.529044450,430270.469574090,429769.576599728,429310.936892201,428863.705546847,428438.619470830,428029.881826754,427632.392707141,427246.184897150,426869.331326223,426512.095248493,423573.372019319,412354.063061564,400116.663216543,387617.450819519,374628.671129686,360544.988590822,345429.709474765,316080.964579299};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={314.509161374653,309.980398017798,309.060421460919,308.725857816145,308.420444716918,308.123579544483,307.842231301249,307.572459413559,307.310835285356,307.057306668409,306.810555274713,306.577130065227,306.371196035435,306.367805075196,306.364074868275,306.360010794064,306.355612294475,306.350868333295,306.345750410671,306.340199085668};

  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=2.0489e6;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=2.0523e6;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={258870.172974984,252026.973076339,247143.028262239,243523.580867570,240944.754557804,238834.783092155,237057.691626821,235406.904532818,233727.820613272,232086.934213579,230504.041464880,228949.748330664,227390.635948594,225890.144270842,224610.308283811,223799.752524797,223833.307373922,224914.147247857,226916.667205483,229509.847375042};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={306.305836749578,305.212616728642,302.480023442784,300.408366650199,298.909663141358,297.670057781293,296.617014160580,295.631677544244,294.622559359966,293.629869759441,292.666318302904,291.714640039147,290.754621536857,289.825718876379,289.029647827961,288.523730457259,288.544672604267,289.218914193266,290.461707168467,292.058148384314};

  parameter Medium_1.MassFlowRate m_flow_init_ss=0;
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
  each nPorts=1) annotation (Placement(transformation(extent={{34,-100},{14,-76}})));

  Modelica.Fluid.Sources.Boundary_pT airsink_ID[Ncell](
  redeclare each package Medium=Medium_2,
  each nPorts=1) annotation (Placement(transformation(extent={{20,-32},{2,-14}})));

  DynamicVCC.Components.Units.Sensors.T_Superheat superheat(
  redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{144,24},{122,46}})));

  DynamicVCC.Components.Units.Sensors.T_Superheat subcooling(
  redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-72,66},{-94,88}})));
/*
  Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    tableName="BC",
    fileName="C:/Jiacheng Ma/BoundaryCondition/JCI/4_30_2021/BC.mat",                                                                                                                                     columns=2:6) annotation (Placement(transformation(extent={{-148,42},{-170,64}})));
  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,                                                                tableName="Mea",
    fileName="C:/Jiacheng Ma/BoundaryCondition/JCI/4_30_2021/Mea.mat",                                                                                                                                 columns=2:10) annotation (Placement(transformation(extent={{-150,70},{-172,92}})));
    */


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
    p_a_start=p_a_start_LL,
    p_b_start=p_b_start_LL,
    h_init=h_init_LL,
    Tt_init=Tt_init_LL)  annotation (Placement(transformation(
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
    p_a_start=p_a_start_VL,
    h_init=h_init_VL,
    Tt_init=Tt_init_VL)       annotation (Placement(transformation(
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
    nPorts=1) annotation (Placement(transformation(extent={{-114,8},{-98,24}})));

   Modelica.Fluid.Sources.Boundary_ph boundary2(
    redeclare package Medium = Medium_1,
    use_p_in=true,
    use_h_in=true,
    nPorts=1) annotation (Placement(transformation(extent={{-116,-24},{-100,-8}})));
*/

  Modelica.Blocks.Sources.Pulse compressorSpeed(
    amplitude=-58.3,
    width=50,
    period=2000,
    offset=58.3,
    startTime=5000)
                   annotation (Placement(transformation(extent={{38,-36},{58,-16}})));
  Modelica.Blocks.Sources.Pulse m_flow_ID(
    amplitude=-1,
    width=50,
    period=2000,
    offset=1,
    startTime=5000)
                   annotation (Placement(transformation(extent={{-172,-12},{-152,8}})));
  Modelica.Blocks.Sources.Pulse m_flow_OD(
    amplitude=-2.1,
    width=50,
    period=2000,
    offset=2.1,
    startTime=5000)
                   annotation (Placement(transformation(extent={{-172,-48},{-152,-28}})));

  Modelica.Blocks.Continuous.FirstOrder firstOrder(
    T=8,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=58.3) annotation (Placement(transformation(extent={{44,8},{64,28}})));
protected
  SI.Mass charge;
  Real COP;
initial equation
  //charge=3;

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

  exv.opening=u[2];
  //compressor.speed=u[1];

/*
  boundary2.p_in=22.602e5;
  boundary2.h_in=253223;
  boundary1.p_in=22.602e5;

  boundary2.p_in=20.524e5;
  boundary2.h_in=261718;
  boundary1.p_in=20.524e5;
  */
  for i in 1:Ncell loop
    airsource_ID[i].m_flow_in=m_flow_air_ID/Ncell;
    airsource_ID[i].T_in=T_ID;
    airsource_ID[i].X_in={1e-4,1-1e-4};
    airsource_OD[i].m_flow_in=m_flow_air_OD/Ncell;
    airsource_OD[i].T_in=T_OD;
    airsource_OD[i].X_in={1e-4,1-1e-4};
  end for;

  m_flow_air_ID=m_flow_ID.y;
  m_flow_air_OD=m_flow_OD.y;
  T_ID=u[5];//294.5;
  T_OD=u[6];//305;

  connect(OutdoorCoil.ports_a2,airsource_OD.ports[1]);
  connect(OutdoorCoil.ports_b2,airsink_OD.ports[1]);
  connect(IndoorCoil.ports_a2,airsource_ID.ports[1]);
  connect(IndoorCoil.ports_b2,airsink_ID.ports[1]);
  connect(subcooling.port, OutdoorCoil.port_b1) annotation (Line(points={{-83,66},{-83,60},{-22,60}},                   color={0,127,255},
      thickness=0.5));
  connect(liquidLine.port_a1, OutdoorCoil.port_b1) annotation (Line(points={{-80,50},{-80,60},{-22,60}}, color={0,127,255}));
  connect(compressor.port_b, OutdoorCoil.port_a1) annotation (Line(points={{70,49},{70,48},{30,48},{30,60},{22,60}}, color={0,127,255}));
  connect(exv.port_b, IndoorCoil.port_a1) annotation (Line(points={{-80,-54},{-80,-60},{-22,-60}}, color={0,127,255}));
  connect(vaporLine.port_b1, compressor.port_a) annotation (Line(points={{112,-14},{112,22},{118,22},{118,49},{112,49}},
                                                                                                                       color={0,127,255}));
  connect(superheat.port, compressor.port_a) annotation (Line(points={{133,24},{134,24},{134,16},{112,16},{112,22},{118,22},{118,49},{112,49}}, color={0,127,255}));
  connect(IndoorCoil.port_b1, vaporLine.port_a1) annotation (Line(points={{22,-60},{112,-60},{112,-38}}, color={0,127,255}));
  connect(liquidLine.port_b1, exv.port_a) annotation (Line(points={{-80,26},{-80,-26}}, color={0,127,255}));
  connect(compressorSpeed.y, firstOrder.u) annotation (Line(points={{59,-26},{64,-26},{64,-4},{36,-4},{36,18},{42,18}}, color={0,0,127}));
  connect(firstOrder.y, compressor.speed) annotation (Line(points={{65,18},{65,34.825},{72.73,34.825}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(extent={{-180,-100},{160,100}})), Icon(coordinateSystem(extent={{-180,-100},{160,100}})),
    experiment(
      StartTime=1,
      StopTime=100000,
      __Dymola_Algorithm="Dassl"));
end TestCycle_OnOff;
