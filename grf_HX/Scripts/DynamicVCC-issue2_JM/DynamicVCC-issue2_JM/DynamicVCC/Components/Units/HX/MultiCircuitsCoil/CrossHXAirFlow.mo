within DynamicVCC.Components.Units.HX.MultiCircuitsCoil;
partial model CrossHXAirFlow "Dry or moist air flow model for cross flow HX"

  import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
  import Modelica.Fluid.Utilities.regStep;
  import Modelica.Constants.T_zero;

  output SI.HeatFlowRate Q_flows "Total heat transfer rate";

  extends DynamicVCC.Interfaces.PartialTwoPort(
    redeclare final package Medium=Medium,
    final enableReverseFlow=false);

  replaceable package Medium=Modelica.Media.Air.MoistAir;

  parameter Integer n=1 "Number of segments in flow direction";

/*************** Connectors ***************/
  DynamicVCC.Interfaces.HeatPorts_a heatPorts[n]
  "Heat transfer with metal wall";

/********************* Fin coil geometry ******************/
  parameter Real eta_fin_overall[n]=ones(n) "Overall fin efficiency";
  parameter SI.Diameter diameters[n]=ones(n) "Tube outter diameter";
  parameter SI.Length lengths[n]=ones(n) "Flow path length";
  parameter SI.Area surfaceAreas[n] "Heat transfer area";

  input SI.Area crossAreas[n](min=1e-10)
    "Effective cross-sectional area (varying in presence of frost)";

/********************* Constants **********************/
  parameter SI.LewisNumber Le=0.854 "Lewis number";
  parameter SI.SpecificEnthalpy delta_h_ig=2836.6e3 "Latent heat of sublimation";
  parameter Real eps=1e-10;

/********************* Initial conditions ******************/
  parameter Medium.MassFlowRate m_flows_init=1 "Initial air flow rate";
  parameter Medium.Temperature T_out_init=Medium.T_default;
  parameter Medium.MassFraction w_out_init=0.01;

/*************** Variables ***************/
  Medium.Temperature T_a "Air inlet temperature";
  Medium.Temperature T_b(start=T_out_init) "Air exit temperatures";
  Medium.MassFraction w_a "Air inlet humidity per unit mass of dry air";
  Medium.MassFraction w_b(start=w_out_init) "Air outlet humidity per unit mass of dry air";
  SI.Velocity v[n](each min=0) "Air flow velocity";
  Medium.BaseProperties mediums[n+1](each T(start=T_out_init)) "Air property at boundary of each segment";

protected
  SI.Temperature TWall[n] "Metal wall surface temperatures";
  Medium.MassFlowRate m_flows[n+1](each start=m_flows_init,each min=0);
  Medium.EnthalpyFlowRate H_flows[n+1];
  Medium.MassFlowRate m_flows_dehumid[n] "Dehumidification water mass flow";

  Medium.SaturationProperties sat[n];
  Real Ntu[n];
  Real Ntu_mass[n];
  SI.HeatFlowRate Q_flows_sen[n] "Sensible heat transfer rates";
  SI.HeatFlowRate Q_flows_lat[n] "Latent heat transfer rate";
  SI.SpecificEnthalpy delta_h_fg[n] "Latent heat of condensation";
  Medium.MassFraction wSat[n] "Saturated humidity ratio";


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

  //w_a={states_a[i].X[Medium.Water]/states_a[i].X[Medium.Air] for i in 1:n};
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

    T_b[i]=T_b_forced[i];
  end for;
  Q_flows=Q_flows_sen+Q_flows_lat;

/*************** Mass balance *********************/
  ports_a.m_flow + ports_b.m_flow=zeros(n);

/*************** momentum balance *********************/
  ports_a.p - ports_b.p = frictionalPressureDrop.dps;

/********************* Boundary Conditions ************************/
  TWall=heatPorts.T;
  heatPorts.Q_flow=Q_flows;
  m_flow[1]=port_a.m_flow;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=225,
      StopTime=2028,
      Tolerance=0.01,
      __Dymola_Algorithm="Dassl"));
end CrossHXAirFlow;
