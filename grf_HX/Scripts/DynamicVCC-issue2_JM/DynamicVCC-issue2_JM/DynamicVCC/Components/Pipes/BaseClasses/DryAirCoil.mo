within DynamicVCC.Components.Pipes.BaseClasses;
partial model DryAirCoil "Dry coil model of air side sensible heat transfer for cross flow HX"

  import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
  import Modelica.Fluid.Utilities.regStep;

  output SI.HeatFlowRate Q_flows[n] "Total heat transfer rate";

  replaceable package Medium=Modelica.Media.Air.MoistAir;

  parameter Integer n=1 "Number of control volumes";

/*************** Connectors ***************/
  DynamicVCC.Interfaces.FluidPorts_a ports_a[n](redeclare each package Medium = Medium,
  each m_flow(min=0))
  "Inlet ports, flow dirction fixed";
  DynamicVCC.Interfaces.FluidPorts_b ports_b[n](redeclare each package Medium=Medium,
  each m_flow(max=0))
  "Outlet ports, flow direction fixed";
  DynamicVCC.Interfaces.HeatPorts_a heatPorts[n]
  "Heat transfer with metal wall";

/********************* Fin coil geometry ******************/
  parameter Real eta_fin_overall=1 "Overall fin efficiency";
  parameter SI.Diameter diameters[n]=ones(n) "Tube outter diameter";
  parameter SI.Length lengths[n]=ones(n) "Flow path length";
  parameter SI.Area surfaceAreas[n] "Heat transfer area";

  input SI.Area crossAreas[n](min=1e-10)
    "Effective cross-sectional area (with frost)";

/********************* Initial conditions ******************/
  parameter Medium.MassFlowRate m_flows_init[n]=ones(n) "Initial air flow rate";
  parameter Medium.Temperature T_out_init[n]=fill(Medium.T_default,n);
  parameter Medium.MassFraction w_out_init[n]=fill(0.01,n);

/*************** Variables ***************/
  Medium.ThermodynamicState states_a[n] "states of inlet ports";
  Medium.ThermodynamicState states_b[n] "states of outlet ports";
  Medium.Temperature T_a[n]=Medium.temperature(states_a);
  output SI.Temperature T_b[n](start=T_out_init) "Air exit temperatures";
  Medium.MassFraction w_a[n] "Air inlet humidity per unit mass of dry air";
  Medium.MassFraction w_b[n](start=w_out_init) "Air outlet humidity per unit mass of dry air";
  SI.Velocity v[n](each min=0) "Air flow velocity";

protected
  SI.Temperature T_s[n] "Surface temperatures";
  SI.Temperature T_b_forced[n](displayUnit="degC",start=T_out_init) "Air exit temperature under forced convection";
  //SI.Temperature T_b_free[n](displayUnit="degC",start=T_out_init) "Under free convection";
  Medium.MassFlowRate m_flows[n](start=m_flows_init,each min=0);
  Medium.SpecificHeatCapacity cp[n]=Medium.specificHeatCapacityCp(states_a);
  Medium.Density rho[n]=Medium.density(states_a);
  Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states_a);
  Medium.SaturationProperties sat[n];


  Real Ntu[n];
  SI.HeatFlowRate Q_flows_sen[n] "Sensible heat transfer rates";
  SI.HeatFlowRate Q_flows_free[n] "Heat transfer rate under free convection";
  SI.HeatFlowRate Q_flows_forced[n] "Heat transfer rate under forced convection";
  parameter SI.MassFlowRate eps=1e-6;

/******************** Heat Transfer *************************/

public
   replaceable model HeatTransfer=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer;
   // Forced convection
   HeatTransfer heatTransfer(
    redeclare final package Medium=Medium,
    final n=n,
    final states=states_a,
    final surfaceAreas=surfaceAreas,
    final dimensions=diameters,
    final lengths=lengths,
    final vs=v);

   // Free convection
   replaceable model FreeConvection=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer;

   FreeConvection freeConvection(
   redeclare final package Medium=Medium,
    final n=n,
    final states=states_a,
    final surfaceAreas=surfaceAreas,
    final dimensions=diameters,
    final lengths=lengths,
    final vs=v);

/******************** Pressure drop *************************/

  replaceable model FrictionalPressureDrop=DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.ConstantFriction
    constrainedby DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.PartialFrictionalPressureDrop;

  FrictionalPressureDrop frictionalPressureDrop(
  redeclare final package Medium=Medium,
  final n=n,
  final states=states_a,
  final surfaceAreas=surfaceAreas,
  final crossAreas=crossAreas,
  final dimensions=diameters,
  final lengths=lengths,
  final vs=v);

equation

/******************* Properties calculation **********************/
  states_a=Medium.setState_phX(ports_a.p,inStream(ports_a.h_outflow),inStream(ports_a.Xi_outflow));
  states_b={Medium.setState_pTX(ports_b[i].p,T_b[i],{max(0,w_b[i]/(1+w_b[i]))}) for i in 1:n};
  w_a={states_a[i].X[Medium.Water]/states_a[i].X[Medium.Air] for i in 1:n};
  for i in 1:n loop
    sat[i].Tsat=T_a[i];
    sat[i].psat=Medium.saturationPressure(T_a[i]);
  end for;

/****************** Heat Transfer ******************/
  for i in 1:n loop
    v[i]=m_flows[i]/(rho[i]*crossAreas[i]);
    Ntu[i]=heatTransfer.alphas[i]*eta_fin_overall*surfaceAreas[i]/((m_flows[i]+eps)*cp[i]);
    Q_flows_sen[i]=m_flows[i]*cp[i]*(1-exp(-Ntu[i]))*(T_s[i]-T_a[i]);
    T_b_forced[i]=T_s[i]-(T_s[i]-T_a[i])*exp(-Ntu[i]);
    //Q_flows_free[i]=m_flows[i]*cp[i]*(T_b_free[i]-T_a[i]);
    //Q_flows_free[i]=freeConvection.alphas[i]*As[i]*(T_s[i]-T_a[i]);
    Q_flows_free[i]=freeConvection.alphas[i]*surfaceAreas[i]*(T_s[i]-T_a[i]);
    Q_flows[i]=regStep(m_flows[i]-0.01,Q_flows_forced[i],Q_flows_free[i],1e-6);
    //Q_flows[i]=Q_flows_forced[i];
    //Q_flows[i]=smooth(0,noEvent(if log10(freeConvection.fluidMotion[i])>1 then Q_flows_free[i] elseif log10(freeConvection.fluidMotion[i])<-2 then Q_flows_forced[i] else (1-log10(freeConvection.fluidMotion[i]))*Q_flows_forced[i]+(log10(freeConvection.fluidMotion[i])+2)*Q_flows_free[i]));
    //T_b[i]=regStep(freeConvection.fluidMotion[i]-0.1,T_b_free[i],T_b_forced[i],1e-4);
    T_b[i]=T_b_forced[i];
  end for;

/*************** Mass balance *********************/
  ports_a.m_flow + ports_b.m_flow=zeros(n);

/*************** momentum balance *********************/
  ports_a.p - ports_b.p = frictionalPressureDrop.dps;

/********************* Boundary Conditions ************************/
  T_s=heatPorts.T;
  heatPorts.Q_flow=Q_flows;
  m_flows=ports_a.m_flow;
  ports_a.Xi_outflow[Medium.Water]=zeros(n); //Reverse flow not allowed
  ports_a.h_outflow=zeros(n); //Reverse flow not allowed
  ports_b.h_outflow=Medium.specificEnthalpy(states_b);
  ports_b.Xi_outflow={{max(0,w_b[i]/(1+w_b[i]))} for i in 1:n};

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=225,
      StopTime=2028,
      Tolerance=0.01,
      __Dymola_Algorithm="Dassl"));
end DryAirCoil;
