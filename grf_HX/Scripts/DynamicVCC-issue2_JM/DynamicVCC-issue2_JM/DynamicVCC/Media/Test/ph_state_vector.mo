within DynamicVCC.Media.Test;
model ph_state_vector
  extends Modelica.Icons.Example;

  replaceable package Medium=DynamicVCC.Media.R410a_NN;
  replaceable package Medium_CP=DynamicVCC.Media.CoolProp.R410a;

  parameter Integer n=5;

  Medium.ThermodynamicState state[n];
  Medium_CP.ThermodynamicState state_CP[n];

  SI.AbsolutePressure p[n];
  SI.SpecificEnthalpy h[n];

  Medium.Density rho[n]=Medium.density(state);

  Medium.SpecificEnthalpy hl[n]=Medium.bubbleEnthalpy(Medium.setSat_p(state.p));

equation
  for i in 1:n loop
    //state[i]=Medium.setState_ph(p[i],h[i]);
    p[i]=(3+0.42*time)*1e5;
    h[i]=(1.1+time*0.035)*1e5;
  end for;

  state_CP=Medium_CP.setState_ph(p,h);
  state=Medium.setState_ph(p,h);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=100,
      Interval=1,
      Tolerance=0.001,
      __Dymola_Algorithm="Dassl"));
end ph_state_vector;
