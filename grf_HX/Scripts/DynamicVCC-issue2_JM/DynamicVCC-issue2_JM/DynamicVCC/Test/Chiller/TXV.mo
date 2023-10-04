within DynamicVCC.Test.Chiller;
model TXV "Thermo-static expansion valve"

  extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialValve;

  Modelica.Blocks.Interfaces.RealInput p_b "Bulb pressure" annotation (Placement(
        transformation(extent={{-8,-66},{24,-34}}), iconTransformation(
        extent={{-11,-11},{11,11}},
        rotation=90,
        origin={-1,-43})));
  Modelica.Blocks.Interfaces.RealInput Q_dot_m "Motor cooling heat";

  import Modelica.Constants.eps;

protected
  parameter SI.AbsolutePressure p_min=80e3;
  parameter SI.Area A_max=250e-6;
  parameter Real k_spring=42e-5 "Spring constant";
  parameter Real coef[2]={0.01764338, -0.24770373} "Coefficients for Area with lift";
  parameter Real lift_max=0.0245092036131709;
  SI.Pressure delta_p "Pressure difference bulb and suction refrigerant";
  Medium.MassFlowRate m_dot_cl "chiller motor cooling line";
  Medium.MassFlowRate m_dot_v "valve flow rate";
  SI.Area Avalve;
  Real lift;
  Medium.Density rho_in=Medium.density(state_a);

equation
  Cd=0.4;
  delta_p=max(eps,(p_b-port_b.p-p_min));
  lift=min(lift_max,k_spring*delta_p/1e3);
  Avalve=max(eps,min(A_max,coef[1]*lift+coef[2]*lift^2));
  m_dot_v=Cd*Avalve*sqrt(rho_in*max(dp,eps));
  m_dot_cl=Cd*8e-6*sqrt(rho_in*max(dp,eps));
  m_flow=m_dot_v+m_dot_cl;
  (port_b.h_outflow-port_a.h_outflow)*m_flow=Q_dot_m;

  port_a.h_outflow=inStream(port_a.h_outflow);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end TXV;
