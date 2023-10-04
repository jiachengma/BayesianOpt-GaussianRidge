within DynamicVCC.Components.Pipes;
model RefFlow1D "Compute source terms and implement flow and volume cells"

  import DynamicVCC.Components.Types.ModelStructure
  "determine number of flow cells (momentum balances) and boundary conditions at port_b";

  //extending connectors at inlet and outlet
  extends DynamicVCC.Interfaces.PartialTwoPort(port_a(m_flow(start=m_flow_init)),port_b(m_flow(start=-m_flow_init)));

  //Volume cells
  extends DynamicVCC.Components.Pipes.BaseClasses.PartialVolumeCell(final n=Ncell, final fluidVolumes={crossAreas[i]*lengths[i] for i in 1:n});

  // Geometry
  parameter SI.Length lengths[n] "Length of refrigerant flow path";
  parameter SI.Area crossAreas[n] "Cross flow area";
  parameter SI.Area surfaceAreas[n] "Surface area (heat transfer area)";
  parameter SI.Diameter diameters[n] "Pipe diameter";

  // Discretization
  parameter Boolean useLumpedPressure=false
  "=true neglect pressure drop across flow path, use a lumped pressure for all volumes";
  parameter Integer Ncell(min=1)=2 "Number of volume cells";
  parameter ModelStructure modelStructure=ModelStructure.av_b;
  final parameter Integer nFlowCell=if useLumpedPressure then nFLumped else nFDistributed;
  final parameter Integer nFDistributed=if modelStructure==ModelStructure.a_v_b       then n+2 elseif modelStructure==ModelStructure.av_vb       then n else n+1
  "Number of Flow cells for distributed volumes";
  final parameter Integer nFLumped=if modelStructure==ModelStructure.a_v_b       then 3 else 2
  "Number of Flow cells under lumped pressure";
  final parameter Integer iLumped=integer(n/2)+1;

  //Initialization
  parameter Medium.MassFlowRate m_flow_init=system.m_flow_init;

  // Flow quantities
  Medium.MassFlowRate m_flows[n+1](each start=m_flow_init, each min=if enableReverseFlow then -Modelica.Constants.inf else 0)
    "Mass flow rate of each volume cell boundary";

  Medium.EnthalpyFlowRate H_flows[n+1]
  "Enthalpy flow rate of each volume cell boundary";

  SI.Velocity v_m[n]={ 0.5*(m_flows[i]+m_flows[i+1])/rho[i]/crossAreas[i] for i in 1:n}
  "Mean velocity of each volume cell";

  // Tube wall heat transfer
  replaceable model HeatTransfer=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer;

  HeatTransfer heatTransfer(
  redeclare final package Medium=Medium,
  final n=n,
  final states=states,
  final surfaceAreas=surfaceAreas,
  final dimensions=diameters,
  final lengths=lengths,
  final vs=v_m);

  DynamicVCC.Interfaces.HeatPorts_a heatPorts[n];

   // Flow model
  replaceable model FlowModel=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.DetailedFlow
  "Momentum balances at staggered grid";

  FlowModel flowModel(
  redeclare final package Medium=Medium,
  final n=nFlowCell,
  final states=statesFM,
  final vs=vsFM,
  final crossAreas=crossAreasFM,
  final dimensions=dimensionsFM,
  final lengths=lengthsFM,
  final enableReverseFlow=enableReverseFlow,
  final m_flow_init=m_flow_init,
  final p_a_start=p_a_start,
  final p_b_start=p_b_start);

protected
  SI.Length lengthsFM[nFlowCell-1] "Length along flow cells";
  SI.Area crossAreasFM[nFlowCell];
  SI.Length dimensionsFM[nFlowCell];
  Medium.ThermodynamicState state_a "State at port_a";
  Medium.ThermodynamicState state_b "State at port_b";
  Medium.ThermodynamicState statesFM[nFlowCell] "States for flow cells";
  SI.Velocity vsFM[nFlowCell] "Flow velocity at boundary of flow cell";

equation

  // Source terms for mass and energy balances
  Q_flows=heatTransfer.Q_flows;
  m_flows_cell=m_flows[1:n]-m_flows[2:n+1];
  H_flows_cell=H_flows[1:n]-H_flows[2:n+1];
/*
  for i in 2:n loop
    H_flows[i]=semiLinear(m_flows[i], Medium.specificEnthalpy(states[i-1]), Medium.specificEnthalpy(states[i]));
  end for;
  H_flows[1]=semiLinear(port_a.m_flow, inStream(port_a.h_outflow), Medium.specificEnthalpy(states[1]));
  H_flows[n+1]=-semiLinear(port_b.m_flow, inStream(port_b.h_outflow), Medium.specificEnthalpy(states[n]));
  */

  // Note h computed from Medium model is based on static quality in two-phase region
  for i in 2:n loop
    H_flows[i]=semiLinear(m_flows[i], h_flow[i-1], h_flow[i]);
  end for;
  H_flows[1]=semiLinear(port_a.m_flow, inStream(port_a.h_outflow), h_flow[1]);
  H_flows[n+1]=-semiLinear(port_b.m_flow, inStream(port_b.h_outflow), h_flow[n]);


  // Boundary conditions
  m_flows[1]=port_a.m_flow;
  m_flows[n+1]=-port_b.m_flow;
  port_a.h_outflow=h_flow[1];
  port_b.h_outflow=h_flow[n];

  state_a=Medium.setState_ph(port_a.p,inStream(port_a.h_outflow));
  state_b=Medium.setState_ph(port_b.p,inStream(port_b.h_outflow));

  // Model structure for evaluating momentum balances
  if useLumpedPressure then
    if modelStructure <> ModelStructure.a_v_b then       //one flow cell
      lengthsFM[1]=sum(lengths);
      crossAreasFM[1:2]={crossAreas[1],crossAreas[n]};
      dimensionsFM[1:2]={diameters[1],diameters[n]};
    else
      lengthsFM[1:2]=fill(sum(lengths)/2,2);
      crossAreasFM[1:3]=fill(sum(crossAreas)/n,3);
      dimensionsFM[1:3]=fill(sum(diameters)/n,3);
    end if;

    if modelStructure <> ModelStructure.av_vb then
      // all pressures are equal
      p[2:n]=fill(p[1],n-1);
    elseif n>2 then
      // need two pressures
      p[2:iLumped-1]=fill(p[1],iLumped-2);
      p[iLumped:n-1]=fill(p[n],n-iLumped);
    end if;

    if modelStructure == ModelStructure.av_vb then
      port_a.p=p[1];
      port_b.p=p[n];
      statesFM[1]=states[1];
      statesFM[2]=states[n];
      flowModel.m_flows[1]=m_flows[iLumped];
      vsFM[1]=v_m[1:iLumped-1]*lengths[1:iLumped-1]/sum(lengths[1:iLumped-1]);
      vsFM[2]=v_m[iLumped:n]*lengths[iLumped:n]/sum(lengths[iLumped:n]);
    elseif modelStructure == ModelStructure.av_b then
      port_a.p=p[1];
      statesFM[1]=states[iLumped];
      statesFM[2]=state_b;
      m_flows[n+1]=flowModel.m_flows[1];
      vsFM[1]=v_m*lengths/sum(lengths);
      vsFM[2]=m_flows[n+1]/Medium.density(state_b)/crossAreas[n];
    elseif modelStructure == ModelStructure.a_v_b then
      statesFM[1]=state_a;
      statesFM[2]=states[iLumped];
      statesFM[3]=state_b;
      m_flows[1]=flowModel.m_flows[1];
      m_flows[n+1]=flowModel.m_flows[2];
      vsFM[1]=m_flows[1]/Medium.density(state_a)/crossAreas[1];
      vsFM[2]=v_m*lengths/sum(lengths);
      vsFM[3]=m_flows[n+1]/Medium.density(state_b)/crossAreas[n];
    end if;

  else
    if modelStructure == ModelStructure.av_vb then
      m_flows[2:n]=flowModel.m_flows;
      lengthsFM=cat(1,{lengths[1]+0.5*lengths[2]},0.5*(lengths[2:n-2]+lengths[3:n-1]),{0.5*lengths[n-1]+lengths[n]});
      crossAreasFM=crossAreas;
      dimensionsFM=diameters;
      vsFM=v_m;
      statesFM[1:n]=states[1:n];
      port_a.p=p[1];
      port_b.p=p[n];
    elseif modelStructure == ModelStructure.av_b then
      m_flows[2:n+1]=flowModel.m_flows;
      lengthsFM=lengths;
      crossAreasFM=cat(1,crossAreas[1:n],{crossAreas[n]});
      dimensionsFM=cat(1,diameters[1:n],{diameters[n]});
      vsFM[1:n]=v_m;
      vsFM[n+1]=m_flows[n+1]/Medium.density(state_b)/crossAreas[n];
      statesFM[1:n]=states;
      statesFM[n+1]=state_b;
      port_a.p=p[1];
    elseif modelStructure == ModelStructure.a_vb then
      m_flows[1:n]=flowModel.m_flows;
      lengthsFM=lengths;
      crossAreasFM=cat(1,{crossAreas[1]},crossAreas[1:n]);
      dimensionsFM=cat(1,{diameters[1]},diameters[1:n]);
      vsFM[2:n+1]=v_m;
      vsFM[1]=m_flows[1]/Medium.density(state_a)/crossAreas[1];
      statesFM[2:n+1]=states;
      statesFM[1]=state_a;
      port_b.p=p[n];
    elseif modelStructure == ModelStructure.a_v_b then
      m_flows[1:n+1]=flowModel.m_flows;
      lengthsFM=cat(1,{0.5*lengths[1]},0.5*(lengths[1:n-1]+lengths[2:n]),{0.5*lengths[n]});
      crossAreasFM=cat(1,{crossAreas[1]},crossAreas[1:n],{crossAreas[n]});
      dimensionsFM=cat(1,{diameters[1]},diameters[1:n],{diameters[n]});
      vsFM[2:n+1]=v_m;
      vsFM[1]=m_flows[1]/Medium.density(state_a)/crossAreas[1];
      vsFM[n+2]=m_flows[n+1]/Medium.density(state_b)/crossAreas[n];
      statesFM[1]=state_a;
      statesFM[2:n+1]=states;
      statesFM[n+2]=state_b;
    else
      assert(false, "Unknown model structure");
    end if;
  end if;

  connect(heatPorts,heatTransfer.heatPorts);

  annotation (Icon(coordinateSystem(preservesurfaceAreaspectRatio=false)), Diagram(coordinateSystem(preservesurfaceAreaspectRatio=false)));
end RefFlow1D;
