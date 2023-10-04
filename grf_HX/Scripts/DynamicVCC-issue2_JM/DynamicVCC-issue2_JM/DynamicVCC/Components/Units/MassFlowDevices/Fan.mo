within DynamicVCC.Components.Units.MassFlowDevices;
model Fan

  import Modelica.Fluid.Utilities.regStep;
  import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
  replaceable package Medium=Modelica.Media.Air.MoistAir;

  Modelica.Blocks.Interfaces.RealInput T_in(unit="K",min=250,start=271);

  Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](each unit="1");

  Modelica.Blocks.Interfaces.RealInput speed "RPM";

  parameter Integer Ncell=1;
  parameter Medium.AbsolutePressure p_atm=Medium.p_default;
  parameter Medium.Density rho_nominal=1.2;
  parameter Medium.AbsolutePressure p_drop_nominal=6;
  parameter Real eta=0.0756 "Fan efficiency";
  parameter SI.Frequency speed_nominal=493;

  // Initial conditions
  parameter SI.VolumeFlowRate V_dot_init=1;

  // Connector
  DynamicVCC.Interfaces.FluidPorts_b ports[Ncell](redeclare each final package Medium=Medium,
  each m_flow(max=0,start=-V_dot_init*rho_nominal/Ncell));

protected
  Medium.Density rho(start=rho_nominal) "Upstream flow density";
  Medium.SpecificEnthalpy h(start=Medium.h_default) "Upstream flow enthalpy";
  SI.Pressure p_drop(min=0,start=6);
  SI.VolumeFlowRate V_dot(start=V_dot_init);
  SI.VolumeFlowRate V_dot_filter(start=V_dot_init);
  SI.VolumeFlowRate V_dot_nominal(start=V_dot_init) "Linear in fan speed";
  SI.VolumeFlowRate V_dot_curve(start=V_dot_init);
  Medium.MassFlowRate m_flows[Ncell](each start=V_dot_init*rho_nominal/Ncell, each min=0);
  SI.Power Pwr;

equation

  ports[2:Ncell].p=fill(ports[1].p,Ncell-1);
  p_drop=sum(ports.p)/Ncell-p_atm;
  rho=Medium.density(Medium.setState_pTX(ports[1].p,T_in,X_in));
  h=Medium.specificEnthalpy(Medium.setState_pTX(p_atm,T_in,X_in));

  sum(m_flows)=rho*V_dot;

  Pwr=V_dot*p_drop/eta;
  V_dot_filter=regStep(speed-490,V_dot_curve,V_dot_nominal,1e-3);
  der(V_dot)=(V_dot_filter-V_dot)/0.01;
  V_dot_nominal=0.0027*speed;

  // Fan curve
  //V_dot = -8.721e-05*p_drop^3 + 0.006327*p_drop^2 -0.1697*p_drop + 2.269;
  //V_dot = max(1e-10,min(2,0.005487*p_drop^3 -0.06925*p_drop^2 -0.3318*p_drop +5.616));

  //V_dot_curve=min(2,max(1e-7,-0.0009489*p_drop^4 + 0.03018*p_drop^3 -0.3424*p_drop^2 + 1.558*p_drop -1.09));
  V_dot_curve=min(2,max(1e-7,0.008056*p_drop^3-0.1756*p_drop^2+0.9753*p_drop-0.2349));

  // Boundary conditions
  for i in 1:Ncell loop
    ports[i].h_outflow=h;
    ports[i].Xi_outflow[Medium.Water]=X_in[Medium.Water];
    ports[i].m_flow=-max(1e-7,m_flows[i]);
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Ellipse(
          extent={{32,0},{48,82}},
          lineColor={28,108,200},
          lineThickness=0.5),
        Ellipse(
          extent={{32,-82},{48,0}},
          lineColor={28,108,200},
          lineThickness=0.5),
        Line(
          points={{-54,40},{12,40}},
          color={28,108,200},
          thickness=0.5,
          arrow={Arrow.None,Arrow.Filled}),
        Line(
          points={{-52,0},{14,0}},
          color={28,108,200},
          thickness=0.5,
          arrow={Arrow.None,Arrow.Filled}),
        Line(
          points={{-52,-40},{14,-40}},
          color={28,108,200},
          thickness=0.5,
          arrow={Arrow.None,Arrow.Filled})}),                    Diagram(coordinateSystem(
          preserveAspectRatio=false)),
    experiment(
      StartTime=615,
      StopTime=10000,
      __Dymola_Algorithm="Dassl"));
end Fan;
