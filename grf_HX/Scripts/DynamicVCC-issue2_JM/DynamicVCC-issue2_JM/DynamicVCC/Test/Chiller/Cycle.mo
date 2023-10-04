within DynamicVCC.Test.Chiller;
model Cycle

  Modelica.Blocks.Interfaces.RealOutput y[6] "Outputs of pressures, water exit temperatures, superheat, subcooling";


  parameter Real u[9]={16553.39114,19778.81283,8857.02706,39925.0108,8334.20534,16793.3871,22665.0566,340.6385,321.782954};
  parameter SI.Mass Mw_c=u[8];
  parameter SI.Mass Mw_e=u[9];



  inner Components.System system(
  p_max=12e5,
  p_min=2e5,
  h_max=5e5,
  h_min=1e5,
  T_max=340,
  T_min=260,
  EnableReverseFlow=EnableReverseFlow) annotation (Placement(transformation(extent={{-162,48},{-126,84}})));

  package Medium_CP=DynamicVCC.Media.CoolProp.R134a;
  //package Medium_NN=DynamicVCC.Media.R134a_NN;

  package Medium_1=Medium_CP;

  package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater (
  cp_const=4186.8,
  d_const=995,
  lambda_const=0.6);

  /************** Initial conditions **************/
  parameter SI.Pressure p_init_c[Ncell]=fill(9.0556e+05,Ncell);
  parameter SI.SpecificEnthalpy h_init_c[Ncell]={425665,422357.093750000,420374.250000000,419181.812500000,418463.156250000,418029.500000000,417767.562500000,414636.218750000,409043.187500000,399053.312500000,381210.031250000,349339.531250000,292414.468750000,250096.875000000,247136.843750000};
  parameter SI.ThermodynamicTemperature Tt_init_c[Ncell]={309.448150634766,309.167724609375,309.000793457031,308.900909423828,308.840911865234,308.804748535156,308.782958984375,308.891448974609,308.885467529297,308.874816894531,308.855804443359,308.821807861328,308.761108398438,306.224029541016,303.171966552734};
  parameter SI.ThermodynamicTemperature Te_init_c[Ncell]={303.050842285156,304.492584228516,306.432006835938,307.517822265625,308.125732421875,308.466064453125,308.656616210938,308.763305664063,308.772247314453,308.787017822266,308.811492919922,308.852111816406,308.919677734375,309.032379150391,309.221282958984};

  //parameter SI.Pressure p_init_e[Ncell]={381953.343750000,381953.343750000,381953.343750000,381953.343750000,381953.343750000,381953.343750000,381953.343750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000};
  parameter SI.Pressure p_init_e[Ncell]=fill(381953.34375,Ncell);
  parameter SI.SpecificEnthalpy h_init_e[Ncell]={248823.453125000,255712.937500000,263863.593750000,273506.250000000,284914.093750000,298410.187500000,314376.812500000,333487.500000000,356096.468750000,382844.187500000,402964.093750000,408547.500000000,410089.375000000,410514.781250000,410632.125000000};
  parameter SI.ThermodynamicTemperature Tt_init_e[Ncell]={280.904602050781,280.938842773438,280.979370117188,281.027343750000,281.084045410156,281.151153564453,281.230560302734,281.281097412109,281.393524169922,281.526519775391,283.943664550781,287.634765625000,288.660552978516,288.943908691406,289.022094726563};
  parameter SI.ThermodynamicTemperature Te_init_e[Ncell]={289.045104980469,289.027313232422,288.962921142578,288.729644775391,287.889068603516,286.771575927734,285.826995849609,285.028564453125,284.361511230469,283.797668457031,283.321044921875,282.918182373047,282.577667236328,282.289825439453,282.046539306641};

  parameter SI.MassFlowRate m_flows_init[Ncell+1]=fill(2.396393,Ncell+1);
  parameter Boolean SteadyState_init=false;
  parameter SI.Pressure dp_init=5e5;
  parameter SI.Power Pwr_init=80e3;
  parameter Real gamma_init=0.9;
  parameter Real Tb_init=284;


  /************** Numerical **************/
  parameter Integer Ncell=15;
  parameter Real C_b=100 "Bulb time constant";

  import DynamicVCC.Components.Types.ModelStructure;
  parameter Boolean EnableReverseFlow=true;
  parameter Boolean useLumpedPressure=true;
  parameter Boolean use_I_flows=true;
  parameter ModelStructure modelStructure=ModelStructure.av_vb;
  import DynamicVCC.Components.Types.DifferentialState;
  parameter DifferentialState differentialState=DifferentialState.pdh;

   /************** Heat transfer **************/
  SI.CoefficientOfHeatTransfer alpha_cond(start=6e4);
  SI.CoefficientOfHeatTransfer alpha_f_c(start=3e3);
  SI.CoefficientOfHeatTransfer alpha_g_c(start=3e3);
  SI.CoefficientOfHeatTransfer alpha_w_c(start=4e4);

  SI.CoefficientOfHeatTransfer alpha_evap(start=6e4);
  SI.CoefficientOfHeatTransfer alpha_f_e(start=3e3);
  //SI.CoefficientOfHeatTransfer alpha_g_e(start=3e3);
  SI.CoefficientOfHeatTransfer alpha_w_e(start=4e4);
/*
  replaceable model Condensation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=5e4,
  final alpha_cst=alpha_cond);

  replaceable model Evaporation =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=5e4,
  final alpha_cst=alpha_evap);

  replaceable model LiquidPhase_c =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=3e3,
  final alpha_cst=alpha_f_c);

  replaceable model VaporPhase_c =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=3e3,
  final alpha_cst=alpha_g_c);

  replaceable model LiquidPhase_c=DynamicVCC.Test.Chiller.HeatTransfer.SinglePhase (
  final C_sf=3);

  replaceable model VaporPhase_c=DynamicVCC.Test.Chiller.HeatTransfer.SinglePhase (
  final C_sf=3);

  replaceable model HeatTransfer_1_c=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=LiquidPhase_c,
  redeclare model VaporZone=VaporPhase_c,
  redeclare model TwoPhaseZone=Condensation);
  */
  replaceable model HeatTransfer_1_c=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
  final alpha_f=u[2],
  final alpha_tp=u[1],
  final alpha_g=u[3]);


/*
  replaceable model LiquidPhase_e =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=3e3,
  final alpha_cst=alpha_f_e);

  replaceable model VaporPhase_e =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=3e3,
  final alpha_cst=alpha_g_e);


  replaceable model HeatTransfer_1_e=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=LiquidPhase_e,
  redeclare model VaporZone=VaporPhase_e,
  redeclare model TwoPhaseZone=Evaporation);
  */

  replaceable model HeatTransfer_1_e=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
  final alpha_f=u[6],
  final alpha_tp=u[5],
  final alpha_g=1e4);


  replaceable model HeatTransfer_2_c = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.ConstHTC (
  final alpha0=alpha_w_c);

  replaceable model HeatTransfer_2_e = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.ConstHTC (
  final alpha0=alpha_w_e);
/*
  replaceable model HeatTransfer_2_c=DynamicVCC.Test.Chiller.HeatTransfer.WaterCond (
  final C_sf=8);

  replaceable model HeatTransfer_2_e=DynamicVCC.Test.Chiller.HeatTransfer.WaterEvap (
  final C_sf=3);
  */

  /************** Friction **************/
  replaceable model Friction_1_c = DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant(f0=0.01);

  replaceable model Friction_1_e = DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant(f0=0.01);

  /************** Components **************/

  replaceable model RefFlow=DynamicVCC.Test.Chiller.RefFlow1D_UniformPressure;

  Components.Units.HX.ShellTubeHX condenser(
  redeclare final package Medium_1=Medium_1,
  redeclare final package Medium_2=Medium_2,
  redeclare final model RefFlow1D=RefFlow,
  redeclare final model HeatTransfer_1=HeatTransfer_1_c,
  redeclare final model HeatTransfer_2=HeatTransfer_2_c,
  redeclare final model Friction_1=Friction_1_c,
  final Ncell=Ncell,
  As_1=23.9328,
  Ac_1=0.0652,
  L_1=2.4384,
  diameter_1=0.01905,
  As_2=19.5231,
  Ac_2=0.0311,
  L_2=2.4384,
  diameter_2=0.01554,
  M_metalWall=Mw_c,
  cp_metalWall=385,
  modelStructure=modelStructure,
  differentialState=differentialState,
  EnableReverseFlow=EnableReverseFlow,
  useLumpedPressure=useLumpedPressure,
  use_I_flows=use_I_flows,
  SteadyState_init=SteadyState_init,
  p_init=p_init_c,
  h_init=h_init_c,
  Tt_init=Tt_init_c,
  Te_init=Te_init_c,
  m_flows_init=m_flows_init) annotation (Placement(transformation(extent={{28,30},{-28,86}})));

  Components.Units.HX.ShellTubeHX evaporator(
  redeclare final package Medium_1=Medium_1,
  redeclare final package Medium_2=Medium_2,
  redeclare final model RefFlow1D=RefFlow,
  redeclare final model HeatTransfer_1=HeatTransfer_1_e,
  redeclare final model HeatTransfer_2=HeatTransfer_2_e,
  redeclare final model Friction_1=Friction_1_e,
  final Ncell=Ncell,
  As_1=22.3716,
  Ac_1=0.07627,
  L_1=2.4384,
  diameter_1=0.0196,
  As_2=18.331,
  Ac_2=0.03,
  L_2=2.4384,
  diameter_2=0.01606,
  M_metalWall=Mw_e,
  cp_metalWall=385,
  modelStructure=modelStructure,
  differentialState=differentialState,
  EnableReverseFlow=EnableReverseFlow,
  useLumpedPressure=useLumpedPressure,
  use_I_flows=use_I_flows,
  SteadyState_init=SteadyState_init,
  p_init=p_init_e,
  h_init=h_init_e,
  Tt_init=Tt_init_e,
  Te_init=Te_init_e,
  m_flows_init=m_flows_init) annotation (Placement(transformation(extent={{-28,-86},{28,-30}})));

  Compressor compressor(
  redeclare package Medium=Medium_1,
  m_flow_init=m_flows_init[1],
  Pwr_init=Pwr_init) annotation (Placement(transformation(
        extent={{31,-31},{-31,31}},
        rotation=-90,
        origin={111,5})));

  Bulb bulb(
  final C=C_b,
  T_init=Tb_init,
  SteadyState_init=SteadyState_init) annotation (Placement(transformation(extent={{-154,-54},{-134,-34}})));

  Modelica.Fluid.Sensors.Temperature temperature_suc(
  redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{86,-92},{62,-68}})));

  Controllernew controller(
  gamma_init=gamma_init) annotation (Placement(transformation(extent={{14,-14},{34,6}})));

  Modelica.Fluid.Sources.MassFlowSource_T sourceCond(nPorts=1,
  redeclare package Medium=Medium_2,
  final use_m_flow_in=true,
  final use_T_in=true) annotation (Placement(transformation(extent={{-58,8},{-38,28}})));

  Modelica.Fluid.Sources.FixedBoundary sinkCond(nPorts=1,
  redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{68,68},{48,88}})));

  Modelica.Fluid.Sources.MassFlowSource_T sourceEvap(nPorts=1,
  redeclare package Medium=Medium_2,
  final use_m_flow_in=true,
  final use_T_in=true) annotation (Placement(transformation(extent={{66,-52},{46,-32}})));

  Modelica.Fluid.Sources.FixedBoundary sinkEvap(nPorts=1,
  redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{-70,-50},{-50,-30}})));

  Modelica.Fluid.Sensors.Temperature temperature_supply(
  redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{-34,-30},{-10,-6}})));

  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:10);
  Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC.mat",columns=2:10);
  Modelica.Blocks.Sources.CombiTimeTable Caldata(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="data_calibration",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/ChillerSimulationData.mat",columns=2:7);



  TXV txv(
  redeclare package Medium=Medium_1,
    dp_nominal=500000,
  m_flow_init=m_flows_init[1],
  m_flow_nominal=2) annotation (Placement(transformation(
        extent={{21,21},{-21,-21}},
        rotation=90,
        origin={-87,-1})));

  Components.Units.Sensors.T_Superheat subcooling(redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-92,64},{-72,84}})));
  Components.Units.Sensors.T_Superheat superheating(redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{128,-56},{148,-36}})));

  SI.Mass charge;
  Real T_sub;

  Modelica.Fluid.Sensors.Temperature temperature_cond(redeclare package Medium = Medium_2) "Condenser water exit temperature" annotation (Placement(transformation(extent={{72,70},{96,94}})));
equation

  charge=condenser.charge+evaporator.charge;
  /*
  alpha_cond=2e4;
  alpha_f_c=1e4;
  alpha_g_c=1e4;
  alpha_w_c=3e4;

  alpha_evap=1e4;
  alpha_f_e=5e3;
  alpha_g_e=5e3;
  alpha_w_e=3e4;
*/

  alpha_cond=u[1];
  alpha_f_c=u[2];
  alpha_g_c=u[3];
  alpha_w_c=u[4];

  alpha_evap=u[5];
  alpha_f_e=u[6];
  alpha_w_e=u[7];

  y[1]=compressor.port_a.p;
  y[2]=compressor.port_b.p;
  y[3]=temperature_supply.T;
  y[4]=temperature_cond.T;
  y[5]=superheating.T;
  y[6]=T_sub;

  controller.u2=BC.y[9];
  sourceCond.m_flow_in=16.9;
  sourceCond.T_in=BC.y[6];
  sourceEvap.m_flow_in=13.6;
  sourceEvap.T_in=BC.y[5];
  /*
  controller.u2=282.55;
  sourceCond.m_flow_in=16.9;
  sourceCond.T_in=302.95;
  sourceEvap.m_flow_in=13.6;
  sourceEvap.T_in=289.05;
  */
  T_sub=-subcooling.T;

  connect(compressor.Q_dot_m,txv.Q_dot_m);

  connect(compressor.port_b, condenser.port_a1) annotation (Line(points={{111,36},{110,36},{110,58},{28,58}}, color={0,127,255}));
  connect(evaporator.port_b1, compressor.port_a) annotation (Line(points={{28,-58},{110,-58},{110,-30},{111,-30},{111,-26}}, color={0,127,255}));
  connect(evaporator.port_b1, temperature_suc.port) annotation (Line(points={{28,-58},{94,-58},{94,-98},{74,-98},{74,-92}}, color={0,127,255}));
  connect(temperature_suc.T, bulb.u) annotation (Line(points={{65.6,-80},{36,-80},{36,-92},{-166,-92},{-166,-44},{-156,-44}}, color={0,0,127}));
  connect(sourceCond.ports[1], condenser.port_a2) annotation (Line(points={{-38,18},{-36,18},{-36,40.64},{-27.44,40.64}}, color={0,127,255}));
  connect(condenser.port_b2, sinkCond.ports[1]) annotation (Line(points={{28,74.8},{38,74.8},{38,78},{48,78}}, color={0,127,255}));
  connect(sourceEvap.ports[1], evaporator.port_a2) annotation (Line(points={{46,-42},{36,-42},{36,-75.36},{27.44,-75.36}}, color={0,127,255}));
  connect(sinkEvap.ports[1], evaporator.port_b2) annotation (Line(points={{-50,-40},{-50,-41.2},{-28,-41.2}}, color={0,127,255}));
  connect(temperature_supply.T, controller.u1) annotation (Line(points={{-13.6,-18},{-2,-18},{-2,2},{12,2}}, color={0,0,127}));
  connect(controller.y, compressor.gamma) annotation (Line(points={{35,-4},{74,-4},{74,5.31},{84.9083,5.31}}, color={0,0,127}));
  connect(temperature_supply.port, evaporator.port_b2) annotation (Line(points={{-22,-30},{-40,-30},{-40,-41.2},{-28,-41.2}}, color={0,127,255}));
  connect(condenser.port_b1, txv.port_a) annotation (Line(points={{-28,58},{-87,58},{-87,20}}, color={0,127,255}));
  connect(txv.port_b, evaporator.port_a1) annotation (Line(points={{-87,-22},{-87,-58},{-28,-58}}, color={0,127,255}));
  connect(bulb.y, txv.p_b) annotation (Line(points={{-133,-44},{-114,-44},{-114,-0.79},{-96.03,-0.79}}, color={0,0,127}));
  connect(subcooling.port, txv.port_a) annotation (Line(points={{-82,64},{-68,64},{-68,58},{-87,58},{-87,20}}, color={0,127,255}));
  connect(superheating.port, compressor.port_a) annotation (Line(points={{138,-56},{110,-56},{110,-30},{111,-30},{111,-26}}, color={0,127,255}));
  connect(condenser.port_b2, temperature_cond.port) annotation (Line(points={{28,74.8},{40,74.8},{40,62},{84,62},{84,70}}, color={0,127,255}));
                     annotation (Placement(transformation(
        extent={{-25,-25},{25,25}},
        rotation=-90,
        origin={-87,1})),
              Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},{180,100}}), graphics={
        Rectangle(
          extent={{-38,40},{32,74}},
          lineColor={0,0,0},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-36,-96},{34,-62}},
          lineColor={0,0,0},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{82,20},{102,20},{122,-20},{62,-20},{82,20}},
          lineColor={0,0,0},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-90,4},{-100,24},{-80,24},{-80,24},{-90,4}},
          lineColor={0,0,0},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-90,4},{-100,-16},{-80,-16},{-80,-16},{-90,4}},
          lineColor={0,0,0},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(
          points={{32,60},{94,60},{94,20}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{34,-80},{96,-80},{96,-20}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{-4,31},{32,31},{32,-19}},
          color={0,0,0},
          thickness=1,
          origin={-59,28},
          rotation=90),
        Line(
          points={{-19,44},{33,44},{33,-20}},
          color={0,0,0},
          thickness=1,
          origin={-57,-36},
          rotation=180)}),                                                                       Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},{180,100}})),
    experiment(
      StartTime=2000,
      StopTime=9000,
      Interval=10,
      __Dymola_Algorithm="Cvode"));
end Cycle;
