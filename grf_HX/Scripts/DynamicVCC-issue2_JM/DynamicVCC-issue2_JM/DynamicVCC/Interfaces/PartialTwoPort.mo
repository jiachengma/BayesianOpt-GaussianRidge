within DynamicVCC.Interfaces;
partial model PartialTwoPort "Inlet and outlet ports of 1-D component"

  import Modelica.Constants;
  outer DynamicVCC.Components.System system;

  replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;

  parameter Boolean enableReverseFlow = system.enableReverseFlow
    "= false restricts flow direction (port_a -> port_b)";

  DynamicVCC.Interfaces.FluidPort_a port_a(redeclare final package Medium = Medium, m_flow(min=if enableReverseFlow then -Constants.inf else 0)) "nominal inlet connector" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  DynamicVCC.Interfaces.FluidPort_b port_b(redeclare final package Medium = Medium, m_flow(max=if enableReverseFlow then +Constants.inf else 0)) "nominal outlet connector" annotation (Placement(transformation(extent={{110,-10},{90,10}}), iconTransformation(extent={{110,-10},{90,10}})));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end PartialTwoPort;
