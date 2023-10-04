within DynamicVCC.Interfaces;
partial model PartialHeatTransfer "Heat transfer model interface"

  replaceable package Medium=Modelica.Media.Interfaces.PartialTwoPhaseMedium;

  parameter Integer n=1 "Number of control volumes";

  input Medium.ThermodynamicState states[n] "Fluid states";

  input SI.Area surfaceAreas[n] "Heat transfer area";

  output SI.HeatFlowRate[n] Q_flows "Heat flow rates";

  // Heat ports
  DynamicVCC.Interfaces.HeatPorts_a heatPorts[n]
  annotation (Placement(transformation(extent={{-10,74},{10,94}}), iconTransformation(extent={{-10,74},{10,94}})));

  // Fluid temperature
  Medium.Temperature T[n]=Medium.temperature(states);
equation
  Q_flows=heatPorts.Q_flow;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end PartialHeatTransfer;
