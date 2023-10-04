within DynamicVCC.Test.Chiller;
model RefFlow1D_UniformPressure "1D flow neglecting pressure drop and eliminate momentum balances (no flow cells)"

  import DynamicVCC.Components.Types.ModelStructure
  "determine number of flow cells (momentum balances) and boundary conditions at port_b";

  //extending connectors at inlet and outlet
  extends DynamicVCC.Interfaces.PartialTwoPort(port_a(m_flow(start=m_flows_init[1])),port_b(m_flow(start=-m_flows_init[Ncell])));

  //Volume cells
  extends DynamicVCC.Components.Pipes.BaseClasses.PartialVolumeCell(final n=Ncell, final V={Ac[i]*L[i] for i in 1:n});

  // Geometry
  parameter SI.Length L[n] "Length of refrigerant flow path";
  parameter SI.Area Ac[n] "Cross flow area";
  parameter SI.Area As[n] "Surface area (heat transfer area)";
  parameter SI.Diameter diameters[n] "Pipe diameter";

  // Discretization
  parameter Integer Ncell(min=1)=2 "Number of volume cells";

  parameter Boolean useLumpedPressure=false
  "=true neglect pressure drop across flow path, use a lumped pressure for all volumes";
  parameter ModelStructure modelStructure=ModelStructure.av_b;
  final parameter Integer nFlowCell=if useLumpedPressure then nFLumped else nFDistributed;
  final parameter Integer nFDistributed=if modelStructure==ModelStructure.a_v_b       then n+2 elseif modelStructure==ModelStructure.av_vb       then n else n+1
  "Number of Flow cells for distributed volumes";
  final parameter Integer nFLumped=if modelStructure==ModelStructure.a_v_b       then 3 else 2
  "Number of Flow cells under lumped pressure";
  final parameter Integer iLumped=integer(n/2)+1;

  //Initialization
  parameter Medium.MassFlowRate m_flows_init[n+1]=ones(n+1);

  // Flow quantities
  Medium.MassFlowRate m_flows[n+1](start=m_flows_init, each min=if EnableReverseFlow then -Modelica.Constants.inf else 0)
    "Mass flow rate of each volume cell boundary";

  Medium.EnthalpyFlowRate H_flows[n+1]
  "Enthalpy flow rate of each volume cell boundary";

  SI.Velocity v_m[n]={ 0.5*(m_flows[i]+m_flows[i+1])/rho[i]/Ac[i] for i in 1:n}
  "Mean velocity of each volume cell";

  // Tube wall heat transfer
  replaceable model HeatTransfer=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.PartialFlowHeatTransfer;

  HeatTransfer heatTransfer(
  redeclare final package Medium=Medium,
  final n=n,
  final states=states,
  final A=As,
  final dimension=diameters,
  final length=L,
  final v=v_m);

  replaceable model FlowCells=DynamicVCC.Components.Pipes.BaseClasses.StaggeredGridMomentum; // Not used

  DynamicVCC.Interfaces.HeatPorts_a heatPorts[n];

protected
  Medium.ThermodynamicState state_a "State at port_a";
  Medium.ThermodynamicState state_b "State at port_b";
equation

  // Source terms for mass and energy balances
  Q_flows=heatTransfer.Q_flows;
  m_flows_cell=m_flows[1:n]-m_flows[2:n+1];
  H_flows_cell=H_flows[1:n]-H_flows[2:n+1];

  for i in 2:n loop
    H_flows[i]=semiLinear(m_flows[i], Medium.specificEnthalpy(states[i-1]), Medium.specificEnthalpy(states[i]));
  end for;
  H_flows[1]=homotopy(semiLinear(port_a.m_flow, inStream(port_a.h_outflow), Medium.specificEnthalpy(states[1])),
  port_a.m_flow*inStream(port_a.h_outflow));
  H_flows[n+1]=homotopy(-semiLinear(port_b.m_flow, inStream(port_b.h_outflow), Medium.specificEnthalpy(states[n])),
  -port_b.m_flow*Medium.specificEnthalpy(states[n]));

  // Boundary conditions
  m_flows[1]=port_a.m_flow;
  m_flows[n+1]=-port_b.m_flow;
  port_a.h_outflow=h[1];
  port_b.h_outflow=h[n];

  state_a=Medium.setState_ph(port_a.p,inStream(port_a.h_outflow));
  state_b=Medium.setState_ph(port_b.p,inStream(port_b.h_outflow));

  // Uniform pressure
  p[2:n]=fill(p[1],n-1);

  port_a.p=p[1];
  port_b.p=p[n];

  connect(heatPorts,heatTransfer.heatPorts);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end RefFlow1D_UniformPressure;
