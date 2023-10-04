within DynamicVCC.Components.Units.HX;
model FinTubeHX "1-D cross-flow fin-and-tube heat exchanger"

  output SI.ThermodynamicTemperature Ta_out_ave(start=Ta_init);

  // extending refrigerant flow models and tube wall models
  extends DynamicVCC.Components.Units.HX.BaseClasses.PartialHX(final C_metalWall=C_FinTube);

  replaceable package Medium_2=Modelica.Media.Air.MoistAir;

  // Tube and fin geometry
  parameter SI.Mass M_tube=1 "Tube wall mass";
  parameter SI.SpecificHeatCapacity cp_tube=1 "Tube material specific heat";
  parameter SI.Area Ac_2 "Air side cross-sectional area";
  parameter Real eta_fin_overall=1 "Overall fin efficiency";
  parameter SI.Area As_2 "Air side heat transfer area";
  parameter SI.Diameter diameter_2;
  parameter SI.Length L_2;
  parameter SI.Mass M_fin=1 "Fin mass";
  parameter SI.SpecificHeatCapacity cp_fin=1 "Fin material specific heat";
  parameter SI.LewisNumber Le=system.Le "Lewis number";
  parameter SI.HeatCapacity C_FinTube=M_tube*cp_tube+M_fin*cp_fin "Heat capacity of tube and fin assuming uniform temperature";

  //Initial conditions
  parameter SI.MassFlowRate m_flow_a_init=1 "Initial air flow rate";
  parameter SI.ThermodynamicTemperature Ta_init=Medium_2.reference_T "Initial air exit temperature";
  parameter SI.MassFraction w_a_init=1e-4 "Initial air exit humidity ratio";

  /*  Air side  */

  replaceable model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer
      constrainedby DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialHeatTransferCoefficient
  "Air side heat transfer model";

  replaceable model FreeConvection =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantHeatTransfer
  "Air side free convection heat transfer model";

  replaceable model FlowModel_2 = DynamicVCC.Components.Pipes.BaseClasses.FrictionalPressureDrop.ConstantFriction
  "Air side frictional pressure drop model";

  // Air flow
  replaceable model AirFlow=DynamicVCC.Components.Pipes.MoistAirCrossFlow;

  AirFlow airFlow(redeclare final package Medium=Medium_2,
  final Ncell=Ncell,
  redeclare final model HeatTransfer=HeatTransfer_2,
  redeclare final model FreeConvection=FreeConvection,
  redeclare final model FrictionalPressureDrop=FlowModel_2,
  final eta_fin_overall=eta_fin_overall,
  final Le=Le,
  final diameters=fill(diameter_2,Ncell),
  final lengths=fill(L_2,Ncell),
  final surfaceAreas=fill(As_2/Ncell,Ncell),
  final crossAreas=fill(Ac_2/Ncell,Ncell),
  final m_flows_init=fill(m_flow_a_init/Ncell,Ncell),
  final T_out_init=fill(Ta_init,Ncell),
  final w_out_init=fill(w_a_init,Ncell));

  // Air flow connectors
  DynamicVCC.Interfaces.FluidPorts_a ports_a2[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
        extent={{-8,-32},{8,32}},
        rotation=90,
        origin={0,-68})));
  DynamicVCC.Interfaces.FluidPorts_b ports_b2[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
        extent={{-8,-32},{8,32}},
        rotation=90,
        origin={0,68})));

equation

  connect(ports_a2,airFlow.ports_a);
  connect(ports_b2,airFlow.ports_b);
  connect(metalWall.heatPorts_b,airFlow.heatPorts);

  Q_flow_2=sum(airFlow.heatPorts.Q_flow);

  Ta_out_ave=sum(airFlow.T_b)/Ncell;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
      Rectangle(
        extent={{-100,-60},{100,60}},
        lineColor={28,108,200},
        lineThickness=0.5),
      Line(
        points={{-100,0},{100,0}},
        color={0,0,0},
        thickness=1),
      Line(
        points={{80,-60},{42,-20},{80,20},{40,60}},
        color={255,0,0},
        thickness=1),
      Line(
        points={{24,-60},{-16,-20},{22,20},{-16,60}},
        color={255,0,0},
        thickness=1),
      Line(
        points={{-34,-60},{-72,-20},{-34,20},{-74,60}},
        color={255,0,0},
        thickness=1),                   Text(
          extent={{-58,80},{56,106}},
          lineColor={0,0,127},
          lineThickness=0.5,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid,
          textString="%name")}),                                 Diagram(coordinateSystem(preserveAspectRatio=false)));
end FinTubeHX;
