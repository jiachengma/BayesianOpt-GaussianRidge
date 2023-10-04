within DynamicVCC.Components.Units.MassFlowDevices;
package BaseClasses
  extends Modelica.Icons.BasesPackage;
  partial model PartialCompressor "Compressor base class"

   extends DynamicVCC.Interfaces.PartialTwoPort(
   final enableReverseFlow=false);//Reverse flow not allowed for compressor

  /************ Initial Conditions *****************/
    parameter Medium.SpecificEnthalpy h_dis_init=Medium.h_default;
    parameter Medium.SpecificEnthalpy h_suc_init=Medium.h_default;
    parameter Medium.AbsolutePressure p_dis_init=Medium.p_default;
    parameter Medium.AbsolutePressure p_suc_init=Medium.p_default;
    parameter Medium.MassFlowRate m_flow_init=system.m_flow_init;
    parameter SI.Power Pwr_init=m_flow_init*(h_dis_init-h_suc_init);

  /*************** Variables *******************/
    parameter SI.Volume Vs=1 "Displacement";
    Medium.MassFlowRate m_flow(start=m_flow_init);
    Medium.AbsolutePressure p_dis(start=p_dis_init) "Discharge pressure";
    Medium.AbsolutePressure p_suc(start=p_suc_init) "Suction pressure";
    Medium.SpecificEnthalpy h_dis(start=h_dis_init) "Discharge enthalpy";
    Medium.SpecificEnthalpy h_suc(start=h_suc_init) "Suction enthalpy";
    SI.Power Pwr(start=Pwr_init) "Power";
    SI.HeatFlowRate Q_loss "Heat loss";
  protected
    Medium.SaturationProperties sat_suc;
    Medium.ThermodynamicState state_suc;
    Medium.SaturationProperties sat_dis;
    Medium.ThermodynamicState state_dis;
    Medium.Temperature T_dis=Medium.temperature(state_dis);
    Medium.Temperature T_suc=Medium.temperature(state_suc);
    //Medium.SpecificEnthalpy h_dis_is(start=h_dis_init) "Isentropic discharge enthalpy";
    Medium.SpecificEntropy s_suc=Medium.specificEntropy(state_suc);

  equation
    sat_dis=Medium.setSat_p(p_dis);
    sat_suc=Medium.setSat_p(p_suc);
    state_dis=Medium.setState_ph(p_dis,h_dis);
    state_suc=Medium.setState_ph(p_suc,h_suc);
    //h_dis_is=Medium.specificEnthalpy(Medium.setState_ps(p_dis,s_suc));

    /* Boundary conditions */
    port_a.p=p_suc;
    port_b.p=p_dis;
    h_dis=port_b.h_outflow;
    h_suc=inStream(port_a.h_outflow);
    port_a.h_outflow=inStream(port_b.h_outflow);//Isenthalpic at flow reversal
    port_a.m_flow=m_flow;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},{100,120}}),
                                                                  graphics={
                             Polygon(
          points={{60,-40},{60,40},{-60,80},{-60,-80},{60,-40}},
          lineColor={0,0,255},
          lineThickness=1,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid,
            rotation=360),                                                  Text(
            extent={{-56,52},{50,158}},
            textColor={0,0,0},
            textString="%name")}),                                 Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-120},{100,120}})));
  end PartialCompressor;

  partial model PartialIsenthalpicValve "Base model of isenthalpic expansion valve"

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialValve;

  equation

    port_a.h_outflow=inStream(port_b.h_outflow);
    port_b.h_outflow=inStream(port_a.h_outflow);

    annotation (Icon(graphics={
          Polygon(
            points={{0,0},{0,0}},
            lineColor={0,0,0},
            fillColor={255,35,38},
            fillPattern=FillPattern.Solid,
            origin={-70,30},
            rotation=90),                                                   Text(
            extent={{-30,42},{28,100}},
            textColor={0,0,0},
            textString="%name")}));
  end PartialIsenthalpicValve;

  partial model PartialValve

    extends DynamicVCC.Interfaces.PartialTwoPort(
    port_a(m_flow(start=m_flow_init),p(start=p_a_init)),port_b(m_flow(start=-m_flow_init),p(start=p_b_init)));

    /******************* Variables ************************/
    parameter SI.Pressure dp_nominal;
    parameter Medium.AbsolutePressure p_a_init=Medium.reference_p;
    parameter Medium.AbsolutePressure p_b_init=Medium.reference_p;
    parameter Medium.MassFlowRate m_flow_init=system.m_flow_init;
    parameter Medium.MassFlowRate m_flow_nominal=system.m_flow_nominal;
    parameter SI.Area Av=1 "Valve fully opened area";
    Real Cd "Discharge coefficient";
    Medium.MassFlowRate m_flow(start=m_flow_init,min=if enableReverseFlow then -Modelica.Constants.inf else 0);

  protected
    Medium.ThermodynamicState state_a;
    Medium.ThermodynamicState state_b;
    SI.Pressure dp(start=p_a_init-p_b_init);

  equation

    state_a=Medium.setState_phX(port_a.p,inStream(port_a.h_outflow));
    state_b=Medium.setState_phX(port_b.p,inStream(port_b.h_outflow));
    dp=port_a.p-port_b.p;

    /* Mass balance */
    port_a.m_flow + port_b.m_flow = 0;
    m_flow=port_a.m_flow;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Polygon(
            points={{0,30},{-40,-50},{40,-50},{0,30}},
            lineColor={0,0,0},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            origin={30,0},
            rotation=90),
          Polygon(
            points={{-40,51},{40,51},{1.77512e-15,-29},{-40,51}},
            lineColor={0,0,0},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            origin={-29,0},
            rotation=90)}),                                        Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialValve;
end BaseClasses;
