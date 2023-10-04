within DynamicVCC.Examples.GreenspeedASHP;
model Accumulator_dfVf "Lumped accumulator model assuming saturated refrigerant states, liquid density and volume as states"

  import DynamicVCC.Components.Types.Dynamics;

  // extending two ports
  extends DynamicVCC.Interfaces.PartialTwoPort;

  parameter Dynamics energyDynamics=system.energyDynamics;
  parameter Dynamics massDynamics=system.massDynamics;

  parameter Medium.AbsolutePressure p_init=Medium.p_default;
  parameter Medium.MassFlowRate m_flow_init=system.m_flow_init;
  parameter SI.Volume V "Volume";
  parameter SI.Volume V_f_init=0.5*V "Initial volume of liquid fluid";

protected
  Medium.ThermodynamicState dewState;
  Medium.ThermodynamicState bubbleState;
  Medium.SaturationProperties sat=Medium.setSat_p(p);
  Medium.AbsolutePressure p(start=p_init);
  Medium.Density rho_f "Bubble density";
  Medium.Density rho_g=Medium.density(dewState) "Dew density";
  Medium.SpecificEnthalpy h_f=Medium.specificEnthalpy(bubbleState);
  Medium.SpecificEnthalpy h_g=Medium.specificEnthalpy(dewState);
  Medium.SpecificInternalEnergy u_f=Medium.specificInternalEnergy(bubbleState);
  Medium.SpecificInternalEnergy u_g=Medium.specificInternalEnergy(dewState);
  Medium.DerEnthalpyByPressure dhdp_f=Medium.dBubbleEnthalpy_dPressure(sat);
  Medium.DerEnthalpyByPressure dhdp_g=Medium.dDewEnthalpy_dPressure(sat);
  Medium.DerDensityByPressure drhodp_f=Medium.dBubbleDensity_dPressure(sat);
  Medium.DerDensityByPressure drhodp_g=Medium.dDewDensity_dPressure(sat);
  SI.Volume V_g(start=V-V_f_init) "Vapor volume";
  SI.Volume V_f(start=V_f_init) "Liquid volume";
  SI.Mass charge;
  SI.Energy energy;

equation
  dewState=Medium.setDewState(sat);
  bubbleState=Medium.setBubbleState(sat);
  charge=V_f*rho_f+V_g*rho_g;
  energy=V_f*rho_f*u_f+V_g*rho_g*u_g;
  /* Mass balance  */
  (V_f+V_g*drhodp_g/drhodp_f)*der(rho_f) + (rho_f-rho_g)*der(V_f) = port_a.m_flow + port_b.m_flow;

  /* Energy balance  */
   (rho_f*u_f-rho_g*u_g)*der(V_f) + (V_f*h_f-V/drhodp_f+V_f*rho_f*dhdp_f/drhodp_f+V_g*rho_g*dhdp_g/drhodp_f+V_g*h_g*drhodp_g/drhodp_f)*der(rho_f)=port_a.m_flow*inStream(port_a.h_outflow) + port_b.m_flow*port_b.h_outflow;

  der(p)=der(rho_f)/drhodp_f;

  V_f+V_g=V;
  port_a.p=p;
  port_b.p=p;
  port_a.h_outflow=h_g; //Reverse flow not allowed
  //port_b.h_outflow=spliceFunction(h_g,inStream(port_a.h_outflow),V_f/V-0.01,1e-7);
  port_b.h_outflow=h_g;

initial equation
   if massDynamics==Dynamics.SteadyState_init then
     der(p)=0;
   elseif massDynamics==Dynamics.Fixed_init then
     p=p_init;
   end if;

   if energyDynamics==Dynamics.SteadyState_init then
     der(V_f)=0;
   elseif energyDynamics==Dynamics.Fixed_init then
     V_f=V_f_init;
   end if;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},{100,120}}),
                                                                graphics={
        Rectangle(
          extent={{-60,60},{60,-20}},
          lineColor={0,0,0},
          lineThickness=1),
        Rectangle(
          extent={{-60,-20},{60,-100}},
          lineColor={0,0,0},
          lineThickness=1,
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),                                 Text(
          extent={{-58,46},{60,164}},
          textColor={0,0,0},
          textString="%name")}),                                 Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},{100,120}})));
end Accumulator_dfVf;