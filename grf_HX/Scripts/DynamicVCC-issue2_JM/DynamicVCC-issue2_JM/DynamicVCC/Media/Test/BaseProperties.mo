within DynamicVCC.Media.Test;
model BaseProperties
  extends Modelica.Icons.Example;

  replaceable package Medium=DynamicVCC.Media.CoolProp.R410a;

  import ExternalMedia.Common.InputChoice;
  //parameter basePropertiesInputChoice=inputChoice.ph;
  parameter Boolean preferredMediumStates=true;

  Medium.BaseProperties medium(
  preferredMediumStates=preferredMediumStates);
equation
  medium.p=10e5+time*1e5;
  medium.h=2e5+2e4*time;

end BaseProperties;
