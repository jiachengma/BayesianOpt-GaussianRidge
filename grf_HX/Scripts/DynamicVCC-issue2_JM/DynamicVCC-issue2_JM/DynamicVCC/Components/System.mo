within DynamicVCC.Components;
model System "Global properties and settings"

  replaceable package Medium=Modelica.Media.Interfaces.PartialMedium;

  parameter SI.AbsolutePressure p_ambient=101325 "Ambient pressure";

  parameter DynamicVCC.Components.Types.Dynamics massDynamics=
                DynamicVCC.Components.Types.Dynamics.DynamicFree_init;

  parameter DynamicVCC.Components.Types.Dynamics energyDynamics=
                DynamicVCC.Components.Types.Dynamics.DynamicFree_init;

  parameter DynamicVCC.Components.Types.Dynamics momentumDynamics=
                DynamicVCC.Components.Types.Dynamics.DynamicFree_init;

  parameter Medium.MassFlowRate m_flow_init=0;

  parameter Medium.MassFlowRate m_flow_small=1e-4;

  parameter Medium.MassFlowRate m_flow_nominal=1;

  parameter SI.AbsolutePressure dp_small(min=0)=1
  "Pressure drop for regularization of zero flow";

  parameter SI.Temperature T_max=500 "Secondary fluid and metal wall maximum temperature";

  parameter SI.Temperature T_min=0 "Secondary fluid and metal wall minimum temperature";

  parameter SI.LewisNumber Le=0.85 "Lewis number";

  parameter SI.Acceleration g=Modelica.Constants.g_n;

  parameter Boolean enableReverseFlow=true;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(extent={{-80,80},{80,-80}}, lineColor={0,0,0}),      Text(
          extent={{-46,-46},{50,50}},
          textColor={28,108,200},
          textString="%name")}),                                 Diagram(coordinateSystem(preserveAspectRatio=false)));
end System;
