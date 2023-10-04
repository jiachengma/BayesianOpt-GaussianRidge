within DynamicVCC.Media.Test;
model PartialDerivatives
  extends Modelica.Icons.Example;

  replaceable package Medium=DynamicVCC.Media.CoolProp.R410a;
  //replaceable package Medium=DynamicVCC.Media.R410A_ph;

  Modelica.SIunits.DerDensityByPressure drhodp_u "partial derivative of density wrt pressure at constant internal energy";
  Modelica.SIunits.DerEnergyByPressure dudp_rho "partial derivative of internal energy wrt pressure at constant density";

  Real dpdrho_u;
  Real dpdu_rho;

  Medium.ThermodynamicState state;
  Medium.SpecificEnthalpy h;
  Medium.Density rho=Medium.density(state);
  Medium.AbsolutePressure p;
  Medium.DerDensityByPressure drhodp_h=Medium.density_derp_h(state)
      "partial derivative of density wrt pressure";
  Medium.DerDensityByEnthalpy drhodh_p=Medium.density_derh_p(state)
      "partial derivative of density wrt enthalpy";

equation
  state = Medium.setState_ph(p,h);
  drhodp_u=(drhodp_h+drhodh_p/rho)/(1+p/rho^2*drhodh_p);
  dudp_rho=-(drhodp_h+drhodh_p/rho)/drhodh_p;
  dpdrho_u=1/drhodp_u;
  dpdu_rho=1/dudp_rho;

  h=(1.3+0.032*time)*1e5;
  p=1e6;

  annotation (experiment(
      StopTime=100,
      __Dymola_NumberOfIntervals=1000,
      __Dymola_Algorithm="Dassl"));
end PartialDerivatives;
