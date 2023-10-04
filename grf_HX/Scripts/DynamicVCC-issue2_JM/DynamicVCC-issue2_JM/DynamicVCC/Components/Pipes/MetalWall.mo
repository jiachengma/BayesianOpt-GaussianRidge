within DynamicVCC.Components.Pipes;
model MetalWall "1D heat conduction of lumped capacity control volumes"

  outer DynamicVCC.Components.System system;

  parameter Integer Ncell=1 "Number of segments";

  parameter SI.HeatCapacity C[Ncell]=fill(1e3,Ncell);
  parameter SI.Temperature T_init[Ncell]=fill(293.15,Ncell);
  parameter DynamicVCC.Components.Types.Dynamics energyDynamics=system.energyDynamics;

  // Variables
  SI.Temperature T[Ncell];
  SI.HeatFlowRate Q_flows_a[Ncell];
  SI.HeatFlowRate Q_flows_b[Ncell];

/*********************Connectors*********************************/

  DynamicVCC.Interfaces.HeatPorts_a heatPorts_a[Ncell] annotation (Placement(transformation(
            extent={{0,-14},{68,6}}),
        iconTransformation(extent={{-8,-56},{6,-42}})));
  DynamicVCC.Interfaces.HeatPorts_b heatPorts_b[Ncell] annotation (Placement(transformation(
            extent={{0,-14},{68,6}}),
        iconTransformation(extent={{-8,42},{6,56}})));
equation

/* Energy balance */
  for i in 1:Ncell loop
    C[i]*der(T[i]) = Q_flows_a[i] - Q_flows_b[i];
  end for;

/* Boundary conditions */
  Q_flows_a = heatPorts_a.Q_flow;
  Q_flows_b = - heatPorts_b.Q_flow;
  T = heatPorts_a.T;
  T = heatPorts_b.T;

initial equation
  if energyDynamics==DynamicVCC.Components.Types.Dynamics.SteadyState_init then
    der(T)=zeros(Ncell);
  elseif energyDynamics==DynamicVCC.Components.Types.Dynamics.Fixed_init then
    T=T_init;
  end if;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                                               Rectangle(
          extent={{-80,-40},{80,40}},
          lineThickness=1,
          lineColor={238,46,47},
          fillColor={0,0,0},
          fillPattern=FillPattern.Forward), Line(
          points={{42,40}},
          color={0,127,255},
          thickness=0.5)}),                                      Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end MetalWall;
