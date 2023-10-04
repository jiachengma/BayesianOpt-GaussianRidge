within DynamicVCC.Media.Test;
model TestMedia

extends Modelica.Icons.Example;

  replaceable package Medium_1=DynamicVCC.Media.CoolProp.R134a;
  replaceable package Medium_2=Modelica.Media.R134a.R134a_ph;


  Medium_1.ThermodynamicState state_1;
  Medium_2.ThermodynamicState state_2;

  SI.AbsolutePressure p;
  SI.SpecificEnthalpy h;

  Medium_1.Density d_1=Medium_1.density(state_1);
  Medium_1.Temperature T_1=Medium_1.temperature(state_1);
  Medium_1.DerDensityByPressure ddph_1=Medium_1.density_derp_h(state_1);
  Medium_1.DerDensityByEnthalpy ddhp_1=Medium_1.density_derh_p(state_1);
  //Medium_1.ThermalConductivity k_1=Medium_1.thermalConductivity(state_1);


  Medium_2.Density d_2=Medium_2.density(state_2);
  Medium_2.Temperature T_2=Medium_2.temperature(state_2);
  Medium_2.DerDensityByPressure ddph_2=Medium_2.density_derp_h(state_2);
  Medium_2.DerDensityByEnthalpy ddhp_2=Medium_2.density_derh_p(state_2);
  //Medium_2.ThermalConductivity k_2=Medium_2.thermalConductivity(state_2);



algorithm
  //p=(2+0.08*time)*1e5;
  p:=8e5;
  h:=(1.1 + time*0.037)*1e5;

  state_1:=Medium_1.setState_phX(p, h);
  state_2:=Medium_2.setState_phX(p, h);



  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=100,
      Interval=0.1,
      __Dymola_Algorithm="Dassl"));
end TestMedia;
