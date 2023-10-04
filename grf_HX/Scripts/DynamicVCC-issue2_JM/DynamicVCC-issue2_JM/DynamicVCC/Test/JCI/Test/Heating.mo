within DynamicVCC.Test.JCI.Test;
package Heating
  extends Modelica.Icons.ExamplesPackage;

  partial model HeatExchangers_heating

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
    pt=0.015875,
    pl=0.020320,
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
    pt=0.021996,
    pl=0.0254,
    t_fin=0.000114,
    Eta_fin_overall=1,
    Ac_e=0.18151322,
    HTA_e=18.19055,
    M_fin=2.820263986,
    HTA_r=0.9327636549,
    Ac_r=0.000355127,
    M_tube=2.904553663);

    replaceable package Medium_1=DynamicVCC.Media.CoolProp.R410a;
    replaceable package Medium_2=Modelica.Media.Air.MoistAir;

    // Initial conditions
    inner DynamicVCC.Components.System system(
    T_max=340,
    T_min=250,
    m_flow_init=0,
    m_flow_nominal=0.051,
    massDynamics=DynamicVCC.Components.Types.Dynamics.DynamicFree_init,
    energyDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
    momentumDynamics=DynamicVCC.Components.Types.Dynamics.SteadyState,
    enableReverseFlow=true);

    parameter Integer Ncell=15;

    parameter Medium_1.AbsolutePressure p_a_start_OD=1544529.375;
    parameter Medium_1.AbsolutePressure p_b_start_OD=1545189;
    parameter Medium_1.SpecificEnthalpy h_init_OD[Ncell]={435574.843750000,428477.718750000,426944.593750000,414900.687500000,400263.500000000,386486.968750000,373628.375000000,361772.375000000,350759.125000000,340430.437500000,330787.843750000,321297.218750000,312257.875000000,304162.687500000,303463.343750000};
    parameter SI.Temperature Tt_init_OD[Ncell]={304.589996337891,297.747924804688,295.971893310547,295.592712402344,295.416229248047,295.310974121094,295.242919921875,295.190704345703,295.148803710938,295.116455078125,295.083770751953,295.080474853516,295.093780517578,295.123565673828,295.104125976563};

    parameter Medium_1.AbsolutePressure p_a_start_ID=1544997.125;
    parameter Medium_1.AbsolutePressure p_b_start_ID=1542945.25;
    parameter Medium_1.SpecificEnthalpy h_init_ID[Ncell]={315522.343750000,288816.562500000,274441.843750000,263236.937500000,250024.843750000,235198.578125000,230914.609375000,232840.046875000,239352.750000000,324321.281250000,425186.875000000,428256.500000000,428553.625000000,428966.812500000,429369.156250000};
    parameter SI.Temperature Tt_init_ID[Ncell]={295.462524414063,295.522064208984,295.511169433594,295.513183593750,295.519409179688,294.133514404297,292.600006103516,294.084350585938,295.456573486328,294.941436767578,295.466156005859,295.507476806641,295.524658203125,295.528686523438,295.528320312500};

    // Numerics
    import DynamicVCC.Components.Types.ModelStructure;
    import DynamicVCC.Components.Types.DifferentialState;
    parameter ModelStructure modelStructure=ModelStructure.av_vb;
    parameter DifferentialState differentialState=DifferentialState.pdh;
    parameter Boolean useLumpedPressure=false;

    parameter Real u_OD[5]={5e3,5e3,2e4,100,17360.3};
    parameter Real u_ID[4]={1e3,3e4,70,5132.6};
    parameter Real u[:]=cat(1,u_OD,u_ID);
    //parameter Real u[9]={1e3,1e3,1e4,150,17360.3,500,5e3,70,5132.6};

    parameter SI.CoefficientOfHeatTransfer alpha_f_OD=u[1];
    parameter SI.CoefficientOfHeatTransfer alpha_g_OD=u[2];
    parameter SI.CoefficientOfHeatTransfer alpha_tp_OD=u[3];
    parameter SI.CoefficientOfHeatTransfer alpha_a_OD=u[4];
    parameter SI.HeatCapacity C_MetalWall_OD=u[5];
  /*
  parameter SI.CoefficientOfHeatTransfer alpha_f_OD=u_OD[1];
  parameter SI.CoefficientOfHeatTransfer alpha_g_OD=u_OD[2];
  parameter SI.CoefficientOfHeatTransfer alpha_tp_OD=u_OD[3];
  parameter SI.CoefficientOfHeatTransfer alpha_a_OD=u_OD[4];
  parameter SI.HeatCapacity C_MetalWall_OD=u_OD[5];
*/

    parameter SI.CoefficientOfHeatTransfer alpha_f_ID=1000;
    parameter SI.CoefficientOfHeatTransfer alpha_g_ID=u[6];
    parameter SI.CoefficientOfHeatTransfer alpha_tp_ID=u[7];
    parameter SI.CoefficientOfHeatTransfer alpha_a_ID=u[8];
    parameter SI.HeatCapacity C_MetalWall_ID=u[9];

    parameter Real f_OD=0.2;
    parameter Real f_ID=0.1;
   /*
  parameter SI.CoefficientOfHeatTransfer alpha_f_ID=1e3;
  parameter SI.CoefficientOfHeatTransfer alpha_g_ID=u_ID[1];
  parameter SI.CoefficientOfHeatTransfer alpha_tp_ID=u_ID[2];
  parameter SI.CoefficientOfHeatTransfer alpha_a_ID=u_ID[3];
  parameter SI.HeatCapacity C_MetalWall_ID=u_ID[4];
*/
    parameter SI.LewisNumber Le=0.854 "Lewis number";

    // Heat transfer and pressure drop
    /*
  replaceable model HeatTransfer_1_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
  final alpha_f=alpha_f_OD,
  final alpha_tp=alpha_tp_OD,
  final alpha_g=alpha_g_OD);
*/

    replaceable model HeatTransfer_1_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.CorrelationPhaseChange (
    redeclare final model TwoPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
    final alpha_nominal=6e3,final m_flow_nominal=0.051,final k=1, final b=0.5, final crossAreas=fill(Geo_OD.Ac_r,Ncell)),
    redeclare final model VaporPhase =
    DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
    final alpha_nominal=1e3,final m_flow_nominal=0.051,final k=1, final b=0.5,final crossAreas=fill(Geo_OD.Ac_r,Ncell)),
    redeclare final model LiquidPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
    final alpha_nominal=1e3,final m_flow_nominal=0.051,final k=1, final b=0.5,final crossAreas=fill(Geo_OD.Ac_r,Ncell)));

    replaceable model HeatTransfer_2_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
    final alpha0=alpha_a_OD);

    replaceable model FreeConvection_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
    final alpha0=5);
  /*
  replaceable model FlowModel_1_OD = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=f_OD);
*/

    replaceable model FlowModel_1_OD=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
    final dp_nominal=5.8343e4/Ncell,
    final k=0.9809,
    final b=0.7651);

    replaceable model Friction_2_OD =
    DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
    f0=0.2);
  /*
  replaceable model HeatTransfer_1_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
      final alpha_f=alpha_f_ID,
      final alpha_tp=alpha_tp_ID,
      final alpha_g=alpha_g_ID);
*/
    replaceable model HeatTransfer_1_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.CorrelationPhaseChange (
      redeclare final model TwoPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
      final alpha_nominal=3e4,final m_flow_nominal=0.051,final k=1, final b=0.5,final crossAreas=fill(Geo_ID.Ac_r,Ncell)),
      redeclare final model VaporPhase =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
      final alpha_nominal=1e3,final m_flow_nominal=0.051,final k=1, final b=0.5,final crossAreas=fill(Geo_ID.Ac_r,Ncell)),
      redeclare final model LiquidPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.NominalHeatTransfer (
      final alpha_nominal=1e3,final m_flow_nominal=0.051,final k=1, final b=0.5,final crossAreas=fill(Geo_ID.Ac_r,Ncell)));

    replaceable model HeatTransfer_2_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
    final alpha0=alpha_a_ID);

    replaceable model FreeConvection_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
    final alpha0=5);

   /*
  replaceable model FlowModel_1_ID=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=f_ID);
*/

    replaceable model FlowModel_1_ID =DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
    final dp_nominal=7.4601e4/Ncell,
    final k=1.0069,
    final b=2.2141);

    replaceable model Friction_2_ID =
    DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
    f0=0.12);

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
    Ac_1=Geo_OD.Ac_r,
    L_1=Geo_OD.L_circuit,
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
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end HeatExchangers_heating;

  model TestCycle_heating
    extends Modelica.Icons.Example;

    Modelica.Blocks.Interfaces.RealOutput charge annotation (Placement(transformation(extent={{160,70},{132,92}}), iconTransformation(extent={{160,70},{132,92}})));
    Modelica.Blocks.Interfaces.RealOutput y[6];
    Modelica.Blocks.Interfaces.RealOutput COP annotation (Placement(transformation(extent={{160,50},{132,72}}), iconTransformation(extent={{160,50},{132,72}})));

    // Components
    extends DynamicVCC.Test.JCI.Test.Heating.HeatExchangers_heating;

    // Initial conditions
    parameter Real speed_init=0;
    parameter Real opening_init=0.4;

    parameter Medium_1.AbsolutePressure p_a_start_LL=1545189;
    parameter Medium_1.AbsolutePressure p_b_start_LL=1544997.5;
    parameter Medium_1.SpecificEnthalpy h_init_LL[Ncell_piping]={428054.593750000,427985.437500000,427851.250000000,427742.937500000};
    parameter SI.Temperature Tt_init_LL[Ncell_piping]={297.077484130859,297.025238037109,296.926483154297,296.846374511719};

    parameter Medium_1.AbsolutePressure p_a_start_VL=1542945.25;
    parameter Medium_1.AbsolutePressure p_b_start_VL=1542850.875;
    parameter Medium_1.SpecificEnthalpy h_init_VL[Ncell_piping]={431775.468750000,431785.500000000,431787.343750000,431787.156250000};
    parameter SI.Temperature Tt_init_VL[Ncell_piping]={299.792022705078,299.798828125000,299.799560546875,299.798828125000};


    DynamicVCC.Test.JCI.Compressor compressor(
    redeclare package Medium=Medium_1,
    Vs=2.9438902e-05,
    m_flow_init=system.m_flow_init,
    m_flow_nominal=0.051,
    speed_nominal=58.3,
    p_dis_init=p_a_start_OD,
    p_suc_init=p_b_start_ID,
    h_dis_init=4.65e5,
    h_suc_init=4.35e5) annotation (Placement(transformation(extent={{112,12},{70,54}})));

    DynamicVCC.Test.JCI.EXV exv(
    redeclare package Medium=Medium_1,
    Av=3.1416e-6,
    p_a_init=p_b_start_OD,
    p_b_init=p_a_start_ID,
    final dp_nominal=10e5) annotation (Placement(transformation(extent={{14,14},{-14,-14}},
          rotation=90,
          origin={-82,30})));

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
    each nPorts=1) annotation (Placement(transformation(extent={{26,-36},{8,-18}})));

    DynamicVCC.Components.Units.Sensors.T_Superheat superheat(
    redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{150,-4},{128,18}})));

    DynamicVCC.Components.Units.Sensors.T_Superheat subcooling(
    redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-38,-50},{-60,-28}})));
  /*
  Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    tableName="BC",
    fileName="C:/Jiacheng Ma/BoundaryCondition/JCI/4_30_2021/BC.mat",                                                                                                                                     columns=2:6) annotation (Placement(transformation(extent={{-124,36},{-146,58}})));
  Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,                                                                tableName="Mea",
    fileName="C:/Jiacheng Ma/BoundaryCondition/JCI/4_30_2021/Mea.mat",                                                                                                                                 columns=2:10) annotation (Placement(transformation(extent={{-124,70},{-146,92}})));
*/
    Modelica.Blocks.Continuous.FirstOrder openingFilter(
      T=5,
      initType=Modelica.Blocks.Types.Init.InitialOutput,
      y_start=opening_init) annotation (Placement(transformation(extent={{-126,20},{-106,40}})));

    // Connect piping between indoor and outdoor units
    replaceable model HeatTransfer_1_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
      alpha0=50);

    replaceable model HeatTransfer_2_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
     final alpha0=5) "Connect piping heat loss to ambient";

    replaceable model FlowModel_1_LL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
    final m_flow_nominal=0.051,
    final dp_nominal=1.3796e4/Ncell_piping,
    final k=0.9796,
    final b=2);

    replaceable model FlowModel_1_VL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
    final m_flow_nominal=0.051,
    final dp_nominal=2.0953e4/Ncell_piping,
    final k=1.0029,
    final b=2);
    /*
  replaceable model FlowModel_1_piping = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
  final lambda0=f_OD);
*/
    parameter Integer Ncell_piping=4;

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
          origin={-80,-28})));

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
      p_b_start=p_b_start_VL,
      h_init=h_init_VL,
      Tt_init=Tt_init_VL)       annotation (Placement(transformation(
          extent={{12,-12},{-12,12}},
          rotation=270,
          origin={46,-34})));

      SI.MassFlowRate m_flow_air_ID;
      SI.MassFlowRate m_flow_air_OD;

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
          origin={52,6})));



          /*
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
        origin={88,-14})));

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
        origin={49,33})));
*/

      /*
  Modelica.Fluid.Sources.MassFlowSource_h boundary(
    redeclare package Medium = Medium_1,
    use_m_flow_in=true,
    use_h_in=true)
              annotation (Placement(transformation(extent={{-140,30},{-120,50}})));
  Modelica.Fluid.Sources.MassFlowSource_h boundary1(
    redeclare package Medium = Medium_1,
    use_m_flow_in=true,
    use_h_in=false)
    annotation (Placement(transformation(extent={{-146,0},{-126,20}})));
    
    */

    Modelica.Blocks.Continuous.FirstOrder filterCompressor(
      T=10,
      initType=Modelica.Blocks.Types.Init.InitialOutput,
      y_start=0)  annotation (Placement(transformation(extent={{100,66},{80,86}})));
    Modelica.Blocks.Continuous.FirstOrder filterOD(
      T=10,
      initType=Modelica.Blocks.Types.Init.InitialOutput,
      y_start=0)   annotation (Placement(transformation(extent={{112,-90},{92,-70}})));
    Modelica.Blocks.Continuous.FirstOrder filterID(
      T=10,
      initType=Modelica.Blocks.Types.Init.InitialOutput,
      y_start=0) annotation (Placement(transformation(extent={{-124,-94},{-144,-74}})));
    Modelica.Blocks.Sources.Step step1(
      height=45,
      offset=0,
      startTime=200)
                    annotation (Placement(transformation(extent={{130,66},{110,86}})));
    Modelica.Blocks.Sources.Step step2(
      height=1,
      offset=0,
      startTime=190)
                    annotation (Placement(transformation(extent={{-82,-94},{-102,-74}})));
    Modelica.Blocks.Sources.Step step3(
      height=3,
      offset=0,
      startTime=190)
                    annotation (Placement(transformation(extent={{148,-90},{128,-70}})));
    Modelica.Blocks.Continuous.PI PI(
      k=-0.0001,
      T=10,
      initType=Modelica.Blocks.Types.Init.InitialOutput,
      y_start=0.4) annotation (Placement(transformation(extent={{-190,-6},{-170,14}})));
    Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(extent={{-200,-50},{-180,-30}})));
    Modelica.Blocks.Sources.Constant const(k=4) annotation (Placement(transformation(extent={{-232,-50},{-212,-30}})));
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax=1, uMin=0.1) annotation (Placement(transformation(extent={{-160,18},{-140,38}})));
  initial equation
    //charge=3.5;

  equation
    //boundary.m_flow_in=0.05;
    //boundary.h_in=2.55e5;
    //boundary1.m_flow_in=-0.05;

    charge=IndoorCoil.charge+OutdoorCoil.charge+liquidLine.charge+vaporLine.charge;
    COP=IndoorCoil.Q_flow_2/(compressor.Pwr+1e-5);
    compressor.T_amb=290;

    y[1]=compressor.port_b.p;
    y[2]=compressor.port_a.p;
    y[3]=IndoorCoil.Ta_out_ave;
    y[4]=OutdoorCoil.Ta_out_ave;
    y[5]=superheat.T;
    y[6]=subcooling.T;

    //speedFilter.u=BC.y[1];
    //openingFilter.u=BC.y[2];
    for i in 1:Ncell loop
      airsource_ID[i].m_flow_in=m_flow_air_ID/Ncell;
      airsource_ID[i].T_in=294;
      airsource_ID[i].X_in={0.0086,1-0.0086};
      airsource_OD[i].m_flow_in=m_flow_air_OD/Ncell;
      airsource_OD[i].T_in=290;
      airsource_OD[i].X_in={0.0086,1-0.0086};
    end for;

    //filterCompressor.u=if time<500 then 58 elseif time<529 then 58-(time-500)*2 else 0;
    //compressor.speed=smooth(0,noEvent(if time<30 then 0 elseif time<59 then (time-30)*2 else 58));
    m_flow_air_ID=filterID.y;//smooth(0,noEvent(if time<50 then 0 elseif time<60 then (time-50)*0.1 else 1));
    m_flow_air_OD=filterOD.y;//smooth(0,noEvent(if time<50 then 0 elseif time<60 then (time-50)*0.21 else 2.1));
    exv.opening=if time<1000 then 0.4 else openingFilter.y;
    reversingValve.opening=0;//smooth(0,noEvent(if time<630 then 1.0 elseif time<632 then 1-(time-630)*0.5 else 0));

    connect(OutdoorCoil.ports_a2,airsource_OD.ports[1]);
    connect(OutdoorCoil.ports_b2,airsink_OD.ports[1]);
    connect(IndoorCoil.ports_a2,airsource_ID.ports[1]);
    connect(IndoorCoil.ports_b2,airsink_ID.ports[1]);
    connect(superheat.port, compressor.port_a) annotation (Line(points={{139,-4},{138,-4},{138,-8},{118,-8},{118,33},{112,33}}, color={0,127,255},
        thickness=0.5));
    connect(reversingValve.port_b2, OutdoorCoil.port_a1) annotation (Line(points={{58.6,-4},{64,-4},{64,60},{22,60}}, color={0,127,255}));
    connect(IndoorCoil.port_b1, vaporLine.port_a1) annotation (Line(points={{22,-60},{46,-60},{46,-46}}, color={0,127,255}));
    connect(vaporLine.port_b1, reversingValve.port_a2) annotation (Line(points={{46,-22},{45.4,-22},{45.4,-4}}, color={0,127,255}));
    connect(compressor.port_b, reversingValve.port_a) annotation (Line(points={{70,33},{70,32},{52,32},{52,16}}, color={0,127,255}));
    connect(reversingValve.port_b, compressor.port_a) annotation (Line(points={{52,-4},{52,-6},{118,-6},{118,33},{112,33}}, color={0,127,255}));
    connect(step1.y, filterCompressor.u) annotation (Line(points={{109,76},{102,76}}, color={0,0,127}));
    connect(step2.y, filterID.u) annotation (Line(points={{-103,-84},{-122,-84}}, color={0,0,127}));
    connect(step3.y, filterOD.u) annotation (Line(points={{127,-80},{114,-80}}, color={0,0,127}));
    connect(filterCompressor.y, compressor.speed) annotation (Line(points={{79,76},{62,76},{62,18.825},{72.73,18.825}}, color={0,0,127}));
    connect(exv.port_b, liquidLine.port_a1) annotation (Line(points={{-82,16},{-80,16},{-80,-16}}, color={0,127,255}));
    connect(liquidLine.port_b1, IndoorCoil.port_a1) annotation (Line(points={{-80,-40},{-80,-60},{-22,-60}}, color={0,127,255}));
    connect(superheat.T, feedback.u2) annotation (Line(points={{131.3,7},{64,7},{64,-6},{-124,-6},{-124,-58},{-190,-58},{-190,-48}}, color={0,0,127}));
    connect(feedback.y, PI.u) annotation (Line(points={{-181,-40},{-174,-40},{-174,-14},{-196,-14},{-196,4},{-192,4}}, color={0,0,127}));
    connect(PI.y, limiter.u) annotation (Line(points={{-169,4},{-169,18},{-170,18},{-170,28},{-162,28}}, color={0,0,127}));
    connect(limiter.y, openingFilter.u) annotation (Line(points={{-139,28},{-134,28},{-134,30},{-128,30}},
                                                                                       color={0,0,127}));
    connect(const.y, feedback.u1) annotation (Line(points={{-211,-40},{-198,-40}}, color={0,0,127}));
    connect(exv.port_a, OutdoorCoil.port_b1) annotation (Line(points={{-82,44},{-82,60},{-22,60}}, color={0,127,255}));
    connect(subcooling.port, IndoorCoil.port_a1) annotation (Line(points={{-49,-50},{-50,-50},{-50,-60},{-22,-60}}, color={0,127,255}));
    annotation (Diagram(coordinateSystem(extent={{-240,-100},{160,100}})), Icon(coordinateSystem(extent={{-240,-100},{160,100}})),
      experiment(
        StartTime=1,
        StopTime=3000,
        __Dymola_Algorithm="Dassl"));
  end TestCycle_heating;
end Heating;
