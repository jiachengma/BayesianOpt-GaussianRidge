within DynamicVCC.Components.Pipes;
model IncompressibleFlow "Incompressible flow using temperature as a state variable"

  import DynamicVCC.Components.Types.Dynamics;

  extends DynamicVCC.Interfaces.PartialTwoPort(
    final enableReverseFlow=false);

  parameter Integer n=2 "Number of control volumes";

  // Geometry
  parameter SI.Length lengths[n] "Length of refrigerant flow path";
  parameter SI.Area crossAreas[n] "Cross flow area";
  parameter SI.Area surfaceAreas[n] "Surface area (heat transfer area)";
  parameter SI.Diameter diameters[n] "Pipe diameter";

  // Initialization
  parameter Medium.Temperature T_init[n]=fill(Medium.T_default,n);
  parameter Types.Dynamics energyDynamics=system.energyDynamics;

  // Variables
  Medium.ThermodynamicState states[n];
  Medium.Temperature T[n](start=T_init);
  Medium.Density rho[n]=Medium.density(states);
  Medium.SpecificHeatCapacity cp[n]=Medium.specificHeatCapacityCp(states);

  // Heat transfer
  DynamicVCC.Interfaces.HeatPorts_b heatPorts[n];

  replaceable model HeatTransfer=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer;

  HeatTransfer heatTransfer(
  redeclare final package Medium=Medium,
  final n=n,
  final states=states,
  final surfaceAreas=surfaceAreas,
  final dimensions=diameters,
  final lengths=lengths,
  final vs=vs);

protected
  SI.Velocity vs[n];
  Medium.MassFlowRate m_flow=port_a.m_flow;
  Medium.AbsolutePressure p=port_a.p;
  Medium.SpecificEnthalpy h[n]=Medium.specificEnthalpy(states);
  SI.HeatCapacity C[n];
  SI.EnthalpyFlowRate H_flow[n+1];

equation

  connect(heatPorts,heatTransfer.heatPorts);

  // Properties
  for i in 1:n loop
    states[i]=Medium.setState_pT(p,T[i]);
    vs[i]=m_flow/rho[i]/crossAreas[i];
    C[i]=rho[i]*lengths[i]*crossAreas[i]*cp[i];
    H_flow[i+1]=m_flow*h[i];
  end for;
  H_flow[1]=m_flow*inStream(port_a.h_outflow);

  // Energy balances
  C.*der(T)=H_flow[1:n]-H_flow[2:n+1]+heatTransfer.Q_flows;

  // Mass balance
  port_a.m_flow + port_b.m_flow = 0;

  // Boundary condition
  port_b.h_outflow=h[n];
  port_a.h_outflow=inStream(port_a.h_outflow);
  port_b.p=p;

initial equation
  if energyDynamics==Dynamics.SteadyState_init then
    der(T)=zeros(n);
  elseif energyDynamics==Dynamics.Fixed_init then
    T=T_init;
  end if;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end IncompressibleFlow;
