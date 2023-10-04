within DynamicVCC.Test.Chiller;
model Compressor "Chiller centrifugal compressor"

  extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialCompressor;

  Modelica.Blocks.Interfaces.RealInput gamma "Inlet guide vane position" annotation (Placement(
        transformation(extent={{-2,-110},{16,-92}}), iconTransformation(
        extent={{-9,-9},{9,9}},
        rotation=90,
        origin={1,-101})));
  Modelica.Blocks.Interfaces.RealOutput Q_dot_m "Motor heat loss";

  // Map coefficients
  constant Real c[:]={1.4354,0.0054572,-0.0030135,-0.017697,0.000005653};
  constant Real b[:]={7.2058,0.8,0.003};
  constant Real a[:]={-0.26524,7.1149,-23.415,0.04173,-0.00089576};
  constant SI.Efficiency eta_em=0.75 "Motor electro-mechanical efficiency";

protected
  constant Real RLA_limit=110;
  SI.MassFlowRate m_flow_max(start=m_flow_init);
  SI.VolumeFlowRate V_dot;
  SI.SpecificEnergy W_p "Ploytropic work";
  SI.Efficiency eta_p(min=0.01,max=1.0) "Polytropic efficiency";
  Real RLA;
  Medium.Density rho_dis=Medium.density(state_dis);
  Medium.Density rho_suc=Medium.density(state_suc);

equation

  V_dot=m_flow/rho_suc;
  W_p=(p_dis/rho_dis-p_suc/rho_suc)*log(p_dis/p_suc)/log((p_dis/rho_dis)/(p_suc/rho_suc));
  Pwr=m_flow*W_p/(eta_p*eta_em);
  Q_loss=Pwr*(1-eta_em);
  h_dis=h_suc+Pwr*eta_em/m_flow;
  m_flow=gamma*m_flow_max;

  /* Mass balance */
  port_a.m_flow + port_b.m_flow = 0;

  /*   Mass flow rate map  */
  m_flow_max=c[1]+c[2]*p_suc/1e3+c[3]*p_dis/1e3+c[4]*SI.Conversions.to_degC(T_suc)+c[5]*p_suc/1e3*p_dis/1e3;

  /*  Polytropic efficiency map  */
  eta_p=a[1]+a[2]*V_dot+a[3]*V_dot^2+a[4]*W_p/1e3+a[5]*(W_p/1e3)^2; // Wp in [kJ/kg] in the map

  /*  RLA map  */
  RLA=b[1]+b[2]*Pwr/1e3+b[3]*(Pwr/1e3)^2;

  Q_dot_m=Q_loss;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Compressor;
