within DynamicVCC.Media.Test;
model ph_state "Compare R410A_ANN with CoolProp"
  extends Modelica.Icons.Example;

  replaceable package Medium_NN=DynamicVCC.Media.R410a_NN;

  //replaceable package Medium_Media=Modelica.Media.R134a.R134a_ph;

  replaceable package Medium_VC=VCLib.Media.Refrigerants.R410a.R410a_IIR_P1_48_T233_473_Horner;

  //Medium_Media.ThermodynamicState state_media;
  Medium_VC.ThermodynamicState state_VC;
  Medium_NN.ThermodynamicState state_NN;

   /*
  Medium_NN.SpecificEntropy s_NN=Medium_NN.specificEntropy(state_NN);
  Medium_CP.SpecificEntropy s_CP=Medium_CP.specificEntropy(state_CP);
  Medium_NN.SpecificEntropy rho_NN=Medium_NN.density(state_NN);
  Medium_CP.SpecificEntropy rho_CP=Medium_CP.density(state_CP);
  Medium_NN.SpecificInternalEnergy u_NN=Medium_NN.specificInternalEnergy(state_NN);
  Medium_CP.SpecificInternalEnergy u_CP=Medium_CP.specificInternalEnergy(state_CP);
*/
  SI.AbsolutePressure p;
  SI.SpecificEnthalpy h;
  Medium_NN.Density d_NN;
  //Medium_Media.Density d_media;
  Medium_VC.Density d_VC;
  /*
  Medium_NN.Temperature T;
  Medium_Media.Temperature T_media;
  Medium_Media.SpecificHeatCapacity cp_media;
  Medium_CP.SpecificHeatCapacity cp_CP;
  Medium_Media.ThermalConductivity k_media;
  Medium_CP.ThermalConductivity k_CP;
  Medium_Media.DynamicViscosity mu_media;
  Medium_CP.DynamicViscosity mu_CP;
*/
equation
  state_NN=Medium_NN.setState_ph(p,h);
  //state_media=Medium_Media.setState_ph(p,h);
  state_VC=Medium_VC.setState_ph(p,h);

  p=(2+0.08*time)*1e5;
  h=(1.1+time*0.034)*1e5;

  d_NN=Medium_NN.density(state_NN);
  //d_media=Medium_Media.density(state_media);
  d_VC=Medium_VC.density(state_VC);

  /*
  cp_media=Medium_Media.specificHeatCapacityCp(state_media);
  cp_CP=Medium_CP.specificHeatCapacityCp(state_CP);
  k_media=Medium_Media.thermalConductivity(state_media);
  k_CP=Medium_CP.thermalConductivity(state_CP);
  mu_media=Medium_Media.dynamicViscosity(state_media);
  mu_CP=Medium_CP.dynamicViscosity(state_CP);
*/
  annotation (experiment(StopTime=100, __Dymola_Algorithm="Dassl"));
end ph_state;
