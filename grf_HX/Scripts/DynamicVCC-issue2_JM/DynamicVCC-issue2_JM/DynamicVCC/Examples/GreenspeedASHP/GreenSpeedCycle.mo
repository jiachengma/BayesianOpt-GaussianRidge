within DynamicVCC.Examples.GreenspeedASHP;
model GreenSpeedCycle
  extends Modelica.Icons.Example;

  extends DynamicVCC.Examples.GreenspeedASHP.HeatExchangers_StaticMom(
    final alpha_f_ID=u[1],
    final alpha_tp_ID=u[2],
    final alpha_g_ID=u[3],
    final alpha_a_ID=u[4],
    final alpha_tp_OD=u[5],
    final alpha_g_OD=u[6],
    final alpha_a_OD=u[7],
    final C_ID=u[8],
    final C_OD=u[9]);

  Modelica.Blocks.Interfaces.RealOutput y[6] "Outputs of pressures, air exit temperatures, subcooling, power";
  Modelica.Blocks.Interfaces.RealOutput y_mea[6];


  // Calibration parameters
  //parameter Real u[10]={5e3, 5e3, 5e3, 50, 5e3, 5e3, 50, 4173.34, 21714.3, 0.0058};
  parameter Real u[10]={15076.25,8641.25,6008.75,60.625,12248.75,13028.75,86.725,2514.43735,16231.4392,0.0056405};

  parameter SI.Volume V_acc=u[10];


  // Piping parameters
  parameter Integer Ncell_piping=5;
  parameter Real lambda0_LL=0.01 "Liquid-line friction factor";
  parameter Real lambda0_VL=0.01 "Vapor-line friction factor";

  // Initial conditions
  parameter SI.Frequency speed_nominal=53;
  parameter Real opening_init=0.16;
  parameter Real alpha_f_init=0.149;

  // Components
  Compressor compressor(
  redeclare final package Medium=Medium_1,
  Vs=2.3765e-5,
  speed_nominal=speed_nominal,
  UA=16.92) annotation (Placement(transformation(
        extent={{-25,-25},{25,25}},
        rotation=90,
        origin={105,-21})));

  DynamicVCC.Components.Units.MassFlowDevices.Accumulator accumulator(
  redeclare package Medium=Medium_1,
  final V=V_acc,
  alpha_f_init=alpha_f_init,
  p_init=p_b_start_OD) annotation (Placement(transformation(extent={{50,-80},{90,-40}})));

  EXV eXV(
  redeclare package Medium=Medium_1,
  Av=3.14e-6,
  final dp_nominal=10e5,
  opening_init=opening_init) annotation (Placement(transformation(
        extent={{-15,-15},{15,15}},
        rotation=-90,
        origin={-79,-41})));

  // Piping heat transfer and pressure drop

  replaceable model HeatTransfer_1_piping =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    final alpha0=5e3);

  replaceable model HeatTransfer_2_piping = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    final alpha0=1);

  replaceable model FlowModel_1_LiquidLine = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=lambda0_LL);

  replaceable model FlowModel_1_VaporLine = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=lambda0_VL);

  Components.Units.HX.Piping liquidLine(
  redeclare package Medium_1=Medium_1,
  redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
  redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
  redeclare final model FlowModel_1=FlowModel_1_LiquidLine,
  final Ncell=Ncell_piping,
  final modelStructure=modelStructure,
  final differentialState=differentialState,
  d_i=0.0110744,
  d_o=0.0127,
  length=10.668,
  p_a_start=p_b_start_ID,
  h_init=fill(h_init_ID[Ncell],Ncell_piping),
  Tt_init=fill(Tt_init_ID[Ncell],Ncell_piping)) annotation (Placement(transformation(
        extent={{14,-14},{-14,14}},
        rotation=90,
        origin={-80,18})));

  Components.Units.HX.Piping vaporLine(
  redeclare package Medium_1=Medium_1,
  redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
  redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
  redeclare final model FlowModel_1=FlowModel_1_VaporLine,
  final Ncell=Ncell_piping,
  final modelStructure=modelStructure,
  final differentialState=differentialState,
  d_i=0.0110744,
  d_o=0.0127,
  length=10.668,
  p_a_start=p_a_start_ID,
  h_init=fill(h_init_ID[1],Ncell_piping),
  Tt_init=fill(Tt_init_ID[1],Ncell_piping)) annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=90,
        origin={106,42})));

  Modelica.Fluid.Sources.MassFlowSource_T fan_ID[Ncell](
  redeclare each final package Medium=Medium_2,
  each use_m_flow_in=true,
  each use_T_in=true,
  each use_X_in=true,
  each nPorts=1) annotation (Placement(transformation(extent={{-36,12},{-16,32}})));

  Modelica.Fluid.Sources.Boundary_pT sink_ID[Ncell](
  redeclare each package Medium=Medium_2,
  each nPorts=1) annotation (Placement(transformation(extent={{-46,90},{-32,104}})));

  Modelica.Fluid.Sources.MassFlowSource_T fan_OD[Ncell](
  redeclare each final package Medium=Medium_2,
  each use_m_flow_in=true,
  each use_T_in=true,
  each use_X_in=true,
  each nPorts=1) annotation (Placement(transformation(extent={{-40,-36},{-20,-16}})));

  Modelica.Fluid.Sources.Boundary_pT sink_OD[Ncell](
  redeclare each package Medium=Medium_2,
  each nPorts=1) annotation (Placement(transformation(extent={{-46,-108},{-32,-94}})));

  SI.Mass charge;

  Modelica.Blocks.Sources.CombiTimeTable BC(
    tableOnFile=true,
    tableName="BC_filter",
    fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/Load-change/BC_filter.mat",
    columns=2:7,
    smoothness=Modelica.Blocks.Types.Smoothness.ModifiedContinuousDerivative)
                                                                "Boundary conditions" annotation (Placement(transformation(extent={{-122,68},{-102,88}})));
  Modelica.Blocks.Sources.CombiTimeTable Mea(
    tableOnFile=true,
    tableName="Mea",
    fileName="C:/Jiacheng Ma/BoundaryCondition/CarrierGreenspeed/Load-change/Mea.mat",
    columns=2:9,
    smoothness=Modelica.Blocks.Types.Smoothness.ModifiedContinuousDerivative)
                                                                "Measurements" annotation (Placement(transformation(extent={{-122,34},{-102,54}})));

  Components.Units.Sensors.T_Superheat subcooling(
  redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-88,70},{-68,90}})));
initial equation
  charge=6.1; // 13 [lb]

equation

  y[1]=compressor.port_b.p;
  y[2]=compressor.port_a.p;
  y[3]=indoorCoil.Ta_out_ave;
  y[4]=outdoorCoil.Ta_out_ave;
  y[5]=subcooling.T;
  y[6]=compressor.Pwr;

  y_mea[1]=Mea.y[1];
  y_mea[2]=Mea.y[2];
  y_mea[3]=Mea.y[3];
  y_mea[4]=Mea.y[4];
  y_mea[5]=Mea.y[6];
  y_mea[6]=Mea.y[7];



  charge = outdoorCoil.charge+indoorCoil.charge+liquidLine.charge+vaporLine.charge+accumulator.charge;
/*
  for i in 1:Ncell loop
    fan_OD[i].m_flow_in=1.81/Ncell;
    fan_OD[i].T_in=293.15;
    fan_OD[i].X_in={1e-4,1-1e-4};
    fan_ID[i].m_flow_in=0.435/Ncell;
    fan_ID[i].T_in=295.4;
    fan_ID[i].X_in={1e-4,1-1e-4};
  end for;

  compressor.speed=53;
  compressor.T_amb=293.15;
  eXV.opening=0.13;
  */
    for i in 1:Ncell loop
    fan_OD[i].m_flow_in=BC.y[4]/Ncell;
    fan_OD[i].T_in=BC.y[6];
    fan_OD[i].X_in={1e-4,1-1e-4};
    fan_ID[i].m_flow_in=BC.y[3]/Ncell;
    fan_ID[i].T_in=BC.y[5];
    fan_ID[i].X_in={1e-4,1-1e-4};
  end for;

  compressor.speed=BC.y[1];
  compressor.T_amb=BC.y[6];
  eXV.opening=BC.y[2];


  connect(fan_OD.ports[1],outdoorCoil.ports_a2);
  connect(outdoorCoil.ports_b2,sink_OD.ports[1]);
  connect(fan_ID.ports[1],indoorCoil.ports_a2);
  connect(indoorCoil.ports_b2,sink_ID.ports[1]);


  connect(indoorCoil.port_b1, liquidLine.port_a1) annotation (Line(
      points={{-28,63},{-28,62},{-80,62},{-80,32}},
      color={0,127,255},
      thickness=0.5));
  connect(liquidLine.port_b1, eXV.port_a) annotation (Line(
      points={{-80,4},{-79,4},{-79,-26}},
      color={0,127,255},
      thickness=0.5));
  connect(eXV.port_b, outdoorCoil.port_a1) annotation (Line(
      points={{-79,-56},{-79,-63},{-30,-63}},
      color={0,127,255},
      thickness=0.5));
  connect(compressor.port_b, vaporLine.port_a1) annotation (Line(
      points={{105,4},{105,14},{106,14},{106,28}},
      color={0,127,255},
      thickness=0.5));
  connect(vaporLine.port_b1, indoorCoil.port_a1) annotation (Line(
      points={{106,56},{106,63},{30,63}},
      color={0,127,255},
      thickness=0.5));
  connect(subcooling.port, liquidLine.port_a1) annotation (Line(points={{-78,70},{-76,70},{-76,62},{-80,62},{-80,32}}, color={0,127,255}));
  connect(outdoorCoil.port_b1, accumulator.port_a) annotation (Line(points={{28,-63},{28,-60},{50,-60}}, color={0,127,255}));
  connect(accumulator.port_b, compressor.port_a) annotation (Line(points={{90,-60},{90,-62},{105,-62},{105,-46}}, color={0,127,255}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-160,-120},{160,120}})),
                                                                 Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,-120},{160,120}})),
    experiment(
      StartTime=6000,
      StopTime=24000,
      Interval=1,
      __Dymola_Algorithm="Dassl"));
end GreenSpeedCycle;
