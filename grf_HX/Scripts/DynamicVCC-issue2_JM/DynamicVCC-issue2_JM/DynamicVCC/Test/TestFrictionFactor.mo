within DynamicVCC.Test;
model TestFrictionFactor "Test friction factor model using piecewise empirical correlations"
  extends Modelica.Icons.Example;

  import Modelica.Constants.pi;
  import DynamicVCC.Utilities.cubicPolynomial;
  import Modelica.Math;

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  parameter SI.Length length=6;
  parameter SI.Diameter diameter=0.0078994;
  parameter SI.Area corssArea=pi/4*diameter^2;

  Medium.MassFlowRate m_flow;
  Medium.AbsolutePressure p;
  Medium.SpecificEnthalpy h;
  Medium.ThermodynamicState state;
  Medium.Density rho=Medium.density(state);
  Medium.DynamicViscosity mu=Medium.dynamicViscosity(state);
  SI.Pressure dp;


equation
  state=Medium.setState_phX(p,h);
  (m_flow)=DynamicVCC.Components.Pipes.BaseClasses.WallFriction.Correlation_SinglePhase.massFlowRate_dp(
  dp,
  rho,
  mu,
  length,
  diameter,
  corssArea,
  1e-4);

  dp=1e3-time*10;

  p=10e5;
  h=2e5;







  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
end TestFrictionFactor;
