within DynamicVCC.Components.Units.MassFlowDevices;
package Compressor "Package of compressor models"
  model Map10Coefficient "Performance map using AHRI 10-coefficient equation"

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialCompressor;
    import Modelica.Math.Vectors;

    // Compressor speed

    Modelica.Blocks.Interfaces.RealInput speed(quantity="Frequency",final unit="Hz",start=speed_init) annotation (Placement(transformation(extent={{54,14},
            {74,34}}),
            iconTransformation(extent={{90,-86},{76,-72}})));

    input SI.Temperature T_amb;

    //SI.Frequency speed(start=speed_init);
    parameter SI.Frequency speed_init=60;

  /******************* Maps *********************/
    constant Real MASS_map[:,:]=[159.291276008480,1.78973432192720,-2.63300401447600,-0.00783168526800000,0.00655027705640000,0.0237112803532000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05;202.524738043980,2.38662488486400,-2.52904562081600,0.000788056896000000,0.00982541558460000,0.0219667500162000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05;208.227524860960,3.90368191664000,-1.88446018172060,0.0166412097874000,-0.0105906570088000,0.0193042983958000,0.000299828196800000,-7.68342117400000e-05,0.000106149940600000,-8.25809736400000e-05;276.874273089700,3.58040601073760,-2.19637876110400,0.0180275432788000,0.0163756926410000,0.0184776903696000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05;353.334649209520,5.23843535200040,-1.45840248837800,0.0419712739094000,0.0254732990930000,0.0136317733376000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05];
    constant Real POWER_map[:,:]=[-279.683678251600,2.10990406957400,13.8321489984790,-0.109808251721000,0.0201255194287000,-0.0794447884860000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;-24.4061357146600,-0.108011374730000,12.4036190197150,-0.109808251721000,0.0301882786604000,-0.0524247383310000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;198.535720623520,-1.76434550257700,11.2779341447780,-0.109808251721000,0.0395242837934000,-0.0273561290030000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;448.101327449020,-3.19728107192900,10.2163985612780,-0.109808251721000,0.0503137980891000,0.00161536680550000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;1122.45682679060,-4.50942876811700,8.66023324954760,-0.109808251721000,0.0782659091039000,0.0766710684376000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000];
    constant SI.Frequency frequency[:]={0,30,45,60,75,117};
    parameter Real Fcor=0.75;

  protected
    Medium.MassFlowRate m_dot_map[size(frequency,1)];
    //Medium.MassFlowRate m_flow_nominal(start=m_flow_init);
    SI.Power Pwr_map[size(frequency,1)];
    Medium.Temperature Te_sat=Medium.saturationTemperature_sat(sat_suc);
    Medium.Temperature Tc_sat=Medium.saturationTemperature_sat(sat_dis);

    //Medium.Density rho_suc_rate "Suction density at rated condition";
    Medium.Density rho_suc=Medium.density(state_suc);
    Real p_ratio(min=0,start=p_dis_init/p_suc_init) "Ratio of discharge and suction pressures";
    Real f_Qloss(start=0.55,min=0.0,max=1.0) "Heat loss coefficient";
    //SI.Frequency f(start=speed_init,min=0);
    //Real UA "Chamber outter wall heat transfer coefficient";
    Real eta_v(start=0.8);

    Medium.SpecificEnthalpy h_dis_is(start=h_dis_init) "Isentropic discharge enthalpy";
    //Medium_CP.ThermodynamicState state_is;

  equation

    //state_is=Medium_CP.setState_ps(p_dis,s_suc);
    h_dis_is=Medium.specificEnthalpy_ps(p_dis,s_suc);

    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;

    m_dot_map = TransientVCC.Utilities.map10coef(
        MASS_map,
        Te_sat,
        Tc_sat) ./ 7937;
    Pwr_map = TransientVCC.Utilities.map10coef(
        POWER_map,
        Te_sat,
        Tc_sat);
    //m_flow=1.0343*((1+Fcor*(rho_suc/rho_suc_rate-1)).*Vectors.interpolate(frequency,m_dot_map,speed)) "Superheat correction";
    //m_flow_nominal=(0.3399/Medium.density(state_suc)/0.0348-0.3927*speed*60/3200+1.0829)*Vectors.interpolate(frequency,m_dot_map,speed);

    eta_v=max(0.6,min(1.0,-0.0308*p_ratio+0.0041*p_ratio^2-0.0873*speed/53+0.9553));
    m_flow=eta_v*Vs*Medium.density(state_suc).*speed;

    //der(m_flow)=(m_flow_nominal-m_flow)/0.1;
    //rho_suc_rate=Medium.density_pT(p_suc,Conv.from_degF(Conv.to_degF(Te_sat)+20));
    Pwr=Vectors.interpolate(frequency,Pwr_map,speed);

    p_ratio=p_dis/p_suc;
    //Q_loss=UA*(T_dis-T_amb);
    //UA=3.596*p_ratio^2 - 30.39*p_ratio + 80.48;
    f_Qloss=min(0.8,max(0.2,-0.6214*m_flow/0.035+9.2704*T_amb/344.494+0.6276*speed/53-7.0899));
    //f_Qloss=0.56;
    Q_loss=Pwr*f_Qloss;
    Pwr=Q_loss+m_flow*(h_dis-h_suc);

  annotation (                     experiment(
        StartTime=22560,
        StopTime=37000,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end Map10Coefficient;

  partial model Efficiencies "Efficiency maps based compressor model"

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialCompressor;

    Modelica.Blocks.Interfaces.RealInput speed(quantity="Frequency",final unit="Hz",start=speed_init) annotation (Placement(transformation(extent={{54,14},
            {74,34}}),
        iconTransformation(extent={{94,-88},{80,-74}})));

    input SI.Temperature T_amb;

    parameter SI.Frequency speed_init=60;
    parameter SI.Frequency speed_nominal=53;

  protected
    Real eta_v(min=0,max=1,start=0.8) "Volumetric efficiency";
    Real eta_is(min=0,max=1,start=0.7) "Isentropic efficiency";
    Medium.Temperature Te_sat=Medium.saturationTemperature_sat(sat_suc);
    Medium.Temperature Tc_sat=Medium.saturationTemperature_sat(sat_dis);
    Medium.Density rho_suc=Medium.density(state_suc) "Suction density";
    Real p_ratio "Ratio of discharge and suction pressures";
    Real f_Qloss(min=0.0,max=1.0) "Heat loss coefficient";
    Medium.SpecificEnthalpy h_dis_is(start=h_dis_init) "Isentropic discharge enthalpy";
  equation

    h_dis_is=Medium.specificEnthalpy_ps(p_dis,s_suc);
    p_ratio=p_dis/p_suc;
    m_flow=homotopy(actual=eta_v*Vs*rho_suc*speed,simplified=m_flow_nominal/speed_nominal*speed);
    eta_is=(h_dis_is-h_suc)*m_flow/Pwr;

    // Efficiency maps
    eta_v=max(0.6,min(1.0,-0.0308*p_ratio+0.0041*p_ratio^2-0.0873*speed/53+0.9553));
    eta_is=max(0.4,min(0.9,0.1312*p_ratio-0.0134*p_ratio^2+0.5623*speed/53-0.0965*p_ratio*speed/53+0.1153));

    Q_loss=Pwr*f_Qloss;
    f_Qloss=min(0.8,max(0.2,-0.6214*m_flow/0.035+9.2704*T_amb/344.494+0.6276*speed/53-7.0899));

    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;

    // Energy balance
    Pwr= m_flow*(h_dis - h_suc) + Q_loss;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-120,-100},{120,120}}),
                                                                  graphics={
                                          Text(
            extent={{-149,-114},{151,-154}},
            textColor={0,0,255},
            textString="%name")}),                                 Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-120,-100},{120,120}})),
      experiment(
        StartTime=615,
        StopTime=11830,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end Efficiencies;

  partial model PartialEfficiencyMap "Partial efficiency maps based compressor model"

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialCompressor;

    Modelica.Blocks.Interfaces.RealInput speed(quantity="Frequency",final unit="Hz") annotation (Placement(transformation(extent={{54,14},
            {74,34}}),
        iconTransformation(extent={{94,-88},{80,-74}})));

    input SI.Temperature T_amb;

    parameter SI.Frequency speed_nominal=53;
    parameter Medium.MassFlowRate m_flow_nominal=system.m_flow_nominal;

  protected
    Medium.Temperature Te_sat=Medium.saturationTemperature_sat(sat_suc);
    Medium.Temperature Tc_sat=Medium.saturationTemperature_sat(sat_dis);
    Medium.Density rho_suc=Medium.density(state_suc) "Suction density";
    Real p_ratio(min=1) "Ratio of discharge and suction pressures";
  equation

    p_ratio=p_dis/p_suc;

    // Mass balance
    port_a.m_flow + port_b.m_flow = 0;

    // Energy balance
    Pwr= m_flow*(h_dis - h_suc) + Q_loss;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialEfficiencyMap;
end Compressor;
