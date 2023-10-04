within DynamicVCC.Components.Units.MassFlowDevices;
package Valve
  model ElectronicExpansionValve "Isenthalpic electronic expansion valve"

    import Modelica.Fluid.Utilities.regRoot2;
    import Modelica.Fluid.Utilities.regStep;

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialIsenthalpicValve;

    Modelica.Blocks.Interfaces.RealInput opening(min=0.0,max=1.0,start=opening_init) "Valve opening" annotation (Placement(transformation(extent={{36,-14},
              {56,6}}),   iconTransformation(extent={{11,-12},{-11,12}},
          rotation=-90,
          origin={1,-46})));
    parameter Real opening_init=0.5;
    parameter SI.Temperature Tcrit=344.494 "Critical temperature";
    parameter Boolean Const_Cd=true;

  protected
    Real p_ratio(start=5);
    SI.TemperatureDifference T_sc_a "subcooling at port_a";
    SI.TemperatureDifference T_sc_b "subcooling at port_a";

  equation
    //Av=4.85e-06;
    //Cd=4.9992e-06*p_ratio^(-0.0546);
    //Cd=4.3843e-06;
    T_sc_a=abs(Medium.saturationTemperature(Medium.pressure(state_a))-Medium.temperature(state_a));
    T_sc_b=abs(Medium.saturationTemperature(Medium.pressure(state_b))-Medium.temperature(state_b));

    //Cd=min(3,max(0.1,0.4875*p_ratio^2-4.114*p_ratio+10.5));
    if Const_Cd then
      Cd=1.8;
      //Cd=regStep(opening-0.8,7,2,1e-5);
    else
      Cd=max(0,exp(1.3113)*(opening)^0.0473*(T_sc_a/Tcrit)^0.2349);
    end if;
    p_ratio=port_a.p/port_b.p;
    m_flow=Av*Cd*opening*regRoot2(2*dp,dp_nominal*1e-9,Medium.density(state_a),Medium.density(state_b));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Text(
            extent={{-36,20},{36,92}},
            textColor={0,0,0},
            textString="%name")}),                                 Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ElectronicExpansionValve;

  model ReversingValve "Reversing valve (four-way valve) to switch flow directions"
    // port_a->compressor discharge, port_b->compressor suction, port_a2->outdoor coil, port_b2->indoor coil
    import DynamicVCC.Components.Types.HPmode;

    extends DynamicVCC.Interfaces.PartialTwoPort;

    Modelica.Blocks.Interfaces.RealInput opening(start=opening_init,max=1,min=0) annotation (Placement(transformation(extent={{-72,58},{-52,78}}),  iconTransformation(extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-52,90})));

    DynamicVCC.Interfaces.FluidPort_a port_a2(redeclare package Medium=Medium)
    "Inlet of the second pair of ports, connect to evaporator" annotation (Placement(transformation(extent={{90,56},{110,76}}),  iconTransformation(extent={{90,56},{110,76}})));

    DynamicVCC.Interfaces.FluidPort_b port_b2(redeclare package Medium=Medium,m_flow(start=-m_flow_init))
    "Outlet of the first pair of ports, connect to condenser" annotation (Placement(transformation(extent={{90,-76},{110,-56}}),  iconTransformation(extent={{90,-76},{110,-56}})));

    parameter HPmode hpMode=HPmode.heating "Operating model only used for initialization";
    parameter SI.Area A_open "Area of each valve when fully opened";
    //parameter SI.Time tau "Time constant of first-order filter";
    parameter Real C_Hd "High temperature stream heat loss coefficient";
    parameter Real C_Hs "Low temperature stream heat loss coefficient";

    // Initial conditions
    parameter SI.Pressure dp_nominal=0.01e5;
    parameter Medium.AbsolutePressure p_dis_init=Medium.reference_p;
    parameter Medium.AbsolutePressure p_suc_init=Medium.reference_p;
    parameter Medium.MassFlowRate m_flow_init=system.m_flow_init;
    parameter Real opening_init=1.0;

    // Components, 4 check valves grid
    DynamicVCC.Components.Units.MassFlowDevices.Valve.CheckValve_Rev checkValve_dis[2](
      redeclare each final package Medium = Medium,
      each final Av=A_open,
      each final dp_nominal=dp_nominal,
      each final C_H=C_Hd,
      each final T_D=T_dis,
      each final T_S=T_suc,
      each final p_a_init=p_dis_init,
      each final p_b_init=p_dis_init,
      final m_flow_init=m_flow_init_dis) " 1: discharge to condenser, 2: discharge to evaporator";

    DynamicVCC.Components.Units.MassFlowDevices.Valve.CheckValve_Rev checkValve_suc[2](
      redeclare each final package Medium = Medium,
      each final Av=2.2*A_open,
      each final dp_nominal=dp_nominal,
      each final C_H=C_Hs,
      each final T_D=T_dis,
      each final T_S=T_suc,
      each final p_a_init=p_suc_init,
      each final p_b_init=p_suc_init,
      final m_flow_init=m_flow_init_suc) " 1: condenser to suction, 2: evaporator to suction";

  protected
    parameter Medium.MassFlowRate m_flow_init_dis[2]= if hpMode==HPmode.heating       then {m_flow_init,0} else {0,m_flow_init};
    parameter Medium.MassFlowRate m_flow_init_suc[2]= if hpMode==HPmode.heating       then {0,m_flow_init} else {m_flow_init,0};
    Medium.Temperature T_dis;
    Medium.Temperature T_suc;

  equation

    T_dis=Medium.temperature(Medium.setState_ph(port_a.p,inStream(port_a.h_outflow)));
    T_suc=Medium.temperature(Medium.setState_ph(port_b.p,port_b.h_outflow));

    checkValve_dis[1].u=opening;
    checkValve_dis[2].u=1.0-opening;
    checkValve_suc[1].u=1.0-opening;
    checkValve_suc[2].u=opening;

    // Boundary conditions
    connect(port_a,checkValve_dis[1].port_a);
    connect(port_a,checkValve_dis[2].port_a);
    connect(port_b,checkValve_suc[1].port_b);
    connect(port_b,checkValve_suc[2].port_b);
    connect(port_a2,checkValve_dis[2].port_b);
    connect(port_a2,checkValve_suc[2].port_a);
    connect(port_b2,checkValve_dis[1].port_b);
    connect(port_b2,checkValve_suc[1].port_a);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Polygon(
            points={{0,0},{0,0}},
            lineColor={0,0,0},
            fillColor={255,35,38},
            fillPattern=FillPattern.Solid,
            origin={-70,30},
            rotation=90),
          Rectangle(
            extent={{-90,-20},{-10,20}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-10,80},{30,-80}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{30,80},{90,52}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{30,14},{90,-14}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{30,-52},{90,-80}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid)}),                      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ReversingValve;

  model CheckValve_Rev "Check valve model used in reverseing valve model. Only one flow direction allowed port_a->port_b"

    import Modelica.Fluid.Utilities.regRoot;

    Modelica.Blocks.Interfaces.RealInput u(max=1.0,min=0.0);

    input Medium.Temperature T_D "Discharge temperature";

    input Medium.Temperature T_S "Suction temperature";

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialValve;

    parameter Real opening_min=0;
    parameter Real C_H "Heat loss coefficient";
  protected
    Real opening(min=opening_min,max=1.0);

  equation
    Cd=1;
    opening=opening_min+(1.0-opening_min)*u;

    m_flow=opening*Cd*Av*sqrt(Medium.density(state_a))*regRoot(dp,1);
    port_b.h_outflow=inStream(port_a.h_outflow)-C_H*(T_D-T_S);
    port_a.h_outflow=inStream(port_a.h_outflow);

    annotation (experiment(
        StartTime=23900,
        StopTime=25170,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end CheckValve_Rev;

  model thermostaticValve
    extends BaseClasses.PartialIsenthalpicValve;
  end thermostaticValve;
end Valve;
