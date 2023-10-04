within DynamicVCC.Examples.GreenspeedASHP;
partial model HeatExchangers_StaticMom

   DynamicVCC.Components.Units.HX.BaseClasses.FinTubeCoil Geo_OD(
    d_o=0.0074,
    d_i=0.0068,
    cp_tube=385,
    rho_tube=8900,
    cp_fin=900,
    rho_fin=2700,
    L_tube=2.7026,
    N_row=2,
    N_prow=44,
    N_circuits=8,
    pf=0.0012,
    pt=0.0216,
    pl=0.0187,
    t_fin=9.9060e-05,
    Eta_fin_overall=1,
    Ac_e=1.556,
    HTA_e=140.64,
    M_fin=18.068);

    DynamicVCC.Components.Units.HX.BaseClasses.FinTubeCoil Geo_ID(
    d_o=0.01,
    d_i=0.0086,
    cp_tube=900,
    rho_tube=2700,
    cp_fin=900,
    rho_fin=2700,
    L_tube=0.4521,
    N_row=3,
    N_prow=28,
    N_circuits=6,
    pf=0.0016,
    pt=0.0254,
    pl=0.0191,
    t_fin=1.0668e-4,
    Eta_fin_overall=1,
    Ac_e=0.183,
    HTA_e=18.82,
    M_fin=2.54);

    // Media
    //replaceable package Medium_1=DynamicVCC.Media.R410a_NN;
    replaceable package Medium_1=DynamicVCC.Media.CoolProp.R410a;
    replaceable package Medium_2=Modelica.Media.Air.MoistAir;

    inner DynamicVCC.Components.System system(
      redeclare package Medium=Medium_1,
      T_max=340,
      T_min=250,
      m_flow_init=m_flow_init,
      m_flow_nominal=0.0373,
      massDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
      energyDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
      momentumDynamics=DynamicVCC.Components.Types.Dynamics.SteadyState,
      enableReverseFlow=true);

    parameter Integer Ncell=15;

    // Initial conditions
    parameter Medium_1.AbsolutePressure p_a_start_ID=3478225;
    parameter Medium_1.AbsolutePressure p_b_start_ID=3.46e6;
    parameter Medium_1.SpecificEnthalpy h_init_ID[Ncell]={430057.593750000,407610,383491,366034.250000000,351035.281250000,337414.843750000,324700.093750000,312634.562500000,301060.937500000,288353.406250000,276881.687500000,267677.562500000,260374.781250000,254624.734375000,250121.687500000};
    parameter SI.Temperature Tt_init_ID[Ncell]={330.776489257813,326.331390380859,326.328460693359,326.326019287109,326.323913574219,326.322174072266,326.320709228516,326.319549560547,326.318634033203,322.648345947266,317.633575439453,313.238616943359,309.553680419922,306.544250488281,304.127410888672};

    parameter Medium_1.AbsolutePressure p_a_start_OD=1.2613e6;
    parameter Medium_1.AbsolutePressure p_b_start_OD=1.2079e6;
    parameter Medium_1.SpecificEnthalpy h_init_OD[Ncell]={246396.031250000,253054.203125000,259872.437500000,266861.531250000,274050.250000000,281472.843750000,289171.125000000,297197.500000000,305620.031250000,314531.125000000,324063.437500000,334422.593750000,345965.187500000,359434.937500000,377347.46875};
    parameter SI.Temperature Tt_init_OD[Ncell]={288.758636474609,288.705413818359,288.664123535156,288.615478515625,288.559265136719,288.495178222656,288.422882080078,288.341949462891,288.251953125000,288.152252197266,288.042205810547,287.920898437500,287.787139892578,287.639282226563,287.392730712891};

    parameter Medium_1.MassFlowRate m_flow_init=0.05;

    parameter Medium_2.MassFlowRate m_flow_a_init_OD=1.7;
    parameter Medium_2.MassFlowRate m_flow_a_init_ID=0.43;
    parameter Medium_2.Temperature Ta_init_OD=287.92;
    parameter Medium_2.Temperature Ta_init_ID=317.78;

    // Numerics
    import DynamicVCC.Components.Types.ModelStructure;
    import DynamicVCC.Components.Types.DifferentialState;
    parameter ModelStructure modelStructure=ModelStructure.av_b;
    parameter DifferentialState differentialState=DifferentialState.pdh;

    parameter SI.HeatCapacity C_OD=21714.3;
    parameter SI.HeatCapacity C_ID=4173.34;

    // Heat transfer, pressure drop, and slip ratio
    //parameter SI.CoefficientOfHeatTransfer alpha_f_OD=5e3;
    parameter SI.CoefficientOfHeatTransfer alpha_g_OD=5e3;
    parameter SI.CoefficientOfHeatTransfer alpha_tp_OD=5e3;
    parameter SI.CoefficientOfHeatTransfer alpha_a_OD=50;

    parameter SI.CoefficientOfHeatTransfer alpha_f_ID=5e3;
    parameter SI.CoefficientOfHeatTransfer alpha_g_ID=5e3;
    parameter SI.CoefficientOfHeatTransfer alpha_tp_ID=5e3;
    parameter SI.CoefficientOfHeatTransfer alpha_a_ID=50;
/*
    parameter Real lambda_tp_OD=0.1;
    parameter Real lambda_g_OD=0.1;
    parameter Real lambda_f_ID=0.1;
    parameter Real lambda_tp_ID=0.1;
    parameter Real lambda_g_ID=0.1;
*/
    replaceable model HeatTransfer_1_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
      final alpha_f=5e3,
      final alpha_tp=alpha_tp_OD,
      final alpha_g=alpha_g_OD);

    replaceable model HeatTransfer_2_OD = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
      final alpha0=alpha_a_OD);

    replaceable model HeatTransfer_1_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
      final alpha_f=alpha_f_ID,
      final alpha_tp=alpha_tp_ID,
      final alpha_g=alpha_g_ID);

    replaceable model HeatTransfer_2_ID = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer (
      final alpha0=alpha_a_ID);
/*
    replaceable model FlowModel_1_OD = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlowPhaseChange (
      final lambda_f=0.1,
      final lambda_tp=lambda_tp_OD,
      final lambda_g=lambda_g_OD);
      */

    replaceable model FlowModel_1_OD = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
      final dp_nominal=57094/Ncell,
      final k=1,
      final b=1);

    replaceable model FlowModel_2_OD = DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.ConstantFriction (
      final f0=0.1);
/*
    replaceable model FlowModel_1_ID = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlowPhaseChange (
      final lambda_f=lambda_f_ID,
      final lambda_tp=lambda_tp_ID,
      final lambda_g=lambda_g_ID);
      */
    replaceable model FlowModel_1_ID = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
      final dp_nominal=19787.57/Ncell,
      final k=1,
      final b=1);

    replaceable model FlowModel_2_ID = DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.ConstantFriction (
      final f0=0.1);

    replaceable model SlipRatio=DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.Zivi;

    // evaporator
    DynamicVCC.Components.Units.HX.FinTubeHX outdoorCoil(
      redeclare final package Medium_1=Medium_1,
      redeclare final package Medium_2=Medium_2,
      redeclare final model HeatTransfer_1=HeatTransfer_1_OD,
      redeclare final model HeatTransfer_2=HeatTransfer_2_OD,
      redeclare final model FlowModel_1=FlowModel_1_OD,
      redeclare final model FlowModel_2=FlowModel_2_OD,
      redeclare final model SlipRatio=SlipRatio,
      final Ncell=Ncell,
      final modelStructure=modelStructure,
      final differentialState=differentialState,
      final C_FinTube=C_OD,
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
      Tt_init=Tt_init_OD,
      m_flow_a_init=m_flow_a_init_OD,
      Ta_init=Ta_init_OD) annotation (Placement(transformation(extent={{-30,-92},{28,-34}})));

    // condenser
    DynamicVCC.Components.Units.HX.FinTubeHX indoorCoil(
    redeclare final package Medium_1 = Medium_1,
    redeclare final package Medium_2 = Medium_2,
    redeclare final model HeatTransfer_1 = HeatTransfer_1_ID,
    redeclare final model HeatTransfer_2 = HeatTransfer_2_ID,
    redeclare final model FlowModel_1 = FlowModel_1_ID,
    redeclare final model FlowModel_2 = FlowModel_2_ID,
    redeclare final model SlipRatio = SlipRatio,
    final Ncell=Ncell,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final C_FinTube=C_ID,
    As_1=Geo_ID.HTA_r,
    Ac_1=Geo_ID.Ac_r,
    L_1=Geo_ID.L_circuit,
    diameter_1=Geo_ID.d_i,
    Ac_2=Geo_ID.Ac_e,
    As_2=Geo_ID.HTA_e,
    diameter_2=Geo_ID.d_o,
    L_2=Geo_ID.L_fin,
    p_a_start=p_a_start_ID,
    p_b_start=p_b_start_ID,
    h_init=h_init_ID,
    Tt_init=Tt_init_ID,
    m_flow_a_init=m_flow_a_init_ID,
    Ta_init=Ta_init_ID) annotation (Placement(transformation(extent={{30,34},{-28,92}})));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{140,100}})),
                                                                 Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{140,100}})));
end HeatExchangers_StaticMom;
