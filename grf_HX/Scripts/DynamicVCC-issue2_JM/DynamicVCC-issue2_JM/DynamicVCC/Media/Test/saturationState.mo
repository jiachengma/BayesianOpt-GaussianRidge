within DynamicVCC.Media.Test;
model saturationState
  extends Modelica.Icons.Example;

  //replaceable package Medium=R410a_NN;
  replaceable package Medium=DynamicVCC.Media.CoolProp.R410a;

  Medium.SaturationProperties sat=Medium.setSat_p(p);

  Medium.AbsolutePressure p;
  Medium.SpecificEnthalpy hl=Medium.bubbleEnthalpy(sat);
  Medium.SpecificEnthalpy hv=Medium.dewEnthalpy(sat);
  Medium.Density rhol=Medium.bubbleDensity(sat);
  Medium.Density rhov=Medium.dewDensity(sat);


equation

 p=(2+time*0.43)*1e5;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end saturationState;
