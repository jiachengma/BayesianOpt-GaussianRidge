within DynamicVCC.Test.JCI;
model Compressor
  import Modelica.Fluid.Utilities.regPow;

  extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialCompressor;

  Modelica.Blocks.Interfaces.RealInput speed(final quantity="Frequency",final unit="Hz") annotation (Placement(transformation(extent={{54,14},
          {74,34}}),
      iconTransformation(extent={{94,-88},{80,-74}})));

  Modelica.Blocks.Interfaces.RealInput T_amb(final quantity="Temperature",final unit="K");

  parameter SI.Frequency speed_nominal=60;
  parameter Medium.MassFlowRate m_flow_nominal=system.m_flow_nominal;

  parameter SI.ThermalConductance UA=8.5969926 "Shell wall thermal conductance";
  parameter Real n=1.3 "Polytropic coefficient";

protected
  Real eta_v(start=0.9,min=0,max=1) "Volumetric efficiency";
  //Real eta_is(min=0,max=1,start=0.7) "Isentropic efficiency";
  Real eta_p(min=0,max=1,start=0.7) "Polytropic efficiency";
  Medium.Temperature Te_sat=Medium.saturationTemperature_sat(sat_suc);
  Medium.Temperature Tc_sat=Medium.saturationTemperature_sat(sat_dis);
  Medium.Density rho_suc=Medium.density(state_suc) "Suction density";
  Real p_ratio(min=1) "Ratio of discharge and suction pressures";
  //Real f_Qloss(min=0.0,max=1.0) "Heat loss coefficient";
  Medium.SpecificEnthalpy h_dis_is(start=h_dis_init) "Isentropic discharge enthalpy";
  SI.Temperature T_w "Average shell wall temperature";
equation

  h_dis_is=Medium.specificEnthalpy_ps(p_dis,s_suc);
  p_ratio=p_dis/p_suc;

  //eta_is=(h_dis_is-h_suc)*m_flow/Pwr;

  // Efficiency maps
  eta_v=0.0018204*speed-0.02571427*p_ratio+0.8846166; // volumetric efficiency
  //eta_is=-0.16854115*p_ratio + 0.07870242*p_ratio^2 + 0.0152457*speed - 0.00525434*speed*p_ratio + 0.49007683;
  eta_p=0.10126318*p_ratio+0.00797325*speed-0.00219957*p_ratio*speed+0.34086291;

  m_flow=homotopy(actual=eta_v*Vs*rho_suc*speed,
                simplified=m_flow_nominal*speed/speed_nominal);

  Pwr=m_flow*(n/(n - 1)*p_suc/rho_suc*(regPow(p_ratio,(n-1)/n) - 1))/eta_p;

  // Heat loss
  T_w=0.767*T_dis+0.233*T_suc;
  Q_loss=UA*(T_w-T_amb);

  // Mass balance
  port_a.m_flow + port_b.m_flow = 0;

  // Energy balance
  Pwr= m_flow*(h_dis - h_suc) + Q_loss;


  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Compressor;
