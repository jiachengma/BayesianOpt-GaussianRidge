within DynamicVCC.Media;
package CoolProp "Call CoolProp via ExternalMedia library"
  extends Modelica.Icons.MaterialPropertiesPackage;

  package R410a
    extends ExternalMedia.Media.CoolPropMedium(
        mediumName="R410A",
        substanceNames={"R410A|debug=1|calc_transport=0|enable_BICUBIC=1"},
        ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
  end R410a;

  package R134a
    extends ExternalMedia.Media.CoolPropMedium(
      mediumName="R134A",
      substanceNames={"R134A|calc_transport=0|enable_BICUBIC=1"},
      ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
  end R134a;
end CoolProp;
