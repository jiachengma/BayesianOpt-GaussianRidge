within DynamicVCC.Components.Units.HX;
model ShellTubeHX "Counter flow shell-and-tube HX"

  extends DynamicVCC.Components.Units.HX.BaseClasses.PartialHX(final C_metalWall=C_tube);

  // Secondary fluid media model
  replaceable package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater;

  // Geometry
  parameter SI.Area Ac_2 "Secondary fluid cross-sectional area";
  parameter SI.Area As_2 "Secondary fluid surface area";
  parameter SI.Diameter diameter_2;
  parameter SI.Length L_2;
  parameter SI.Mass M_metalWall "Tube wall mass";
  parameter SI.SpecificHeatCapacity cp_metalWall "Tube material specific heat";
  parameter SI.HeatCapacity C_tube=M_metalWall*cp_metalWall;

  // Initial conditions
  parameter Medium_2.Temperature Te_init[Ncell]=fill(Medium_2.T_default,Ncell)
   "Initial conditions of the secondary fluid temperature";

  // Secondary fluid heat transfer
  replaceable model HeatTransfer_2 = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer;

  /* Secondary fluid */
  replaceable model SecondaryFluid =
      DynamicVCC.Components.Pipes.IncompressibleFlow;

  SecondaryFluid secondaryFluid(
  redeclare final package Medium=Medium_2,
  redeclare final model HeatTransfer=HeatTransfer_2,
  final n=Ncell,
  final lengths=fill(L_2/Ncell,Ncell),
  final crossAreas=fill(Ac_2,Ncell),
  final surfaceAreas=fill(As_2/Ncell,Ncell),
  final diameters=fill(diameter_2,Ncell),
  final T_init=Te_init);

  // Secondary fluid connectors
  DynamicVCC.Interfaces.FluidPort_a port_a2(redeclare final package Medium=Medium_2) annotation (Placement(transformation(extent={{88,-72},{108,-52}}), iconTransformation(extent={{88,-72},{108,-52}})));
  DynamicVCC.Interfaces.FluidPort_b port_b2(redeclare final package Medium=Medium_2) annotation (Placement(transformation(extent={{-110,50},{-90,70}}), iconTransformation(extent={{-110,50},{-90,70}})));

equation
  connect(port_a2,secondaryFluid.port_a);
  connect(port_b2,secondaryFluid.port_b);
  for i in 1:Ncell loop
    connect(secondaryFluid.heatPorts[i],metalWall.heatPorts_b[Ncell-i+1]);
  end for;

  Q_flow_2=sum(secondaryFluid.heatTransfer.Q_flows);

  annotation (Icon(graphics={
        Rectangle(
          extent={{-100,30},{100,-30}},
          pattern=LinePattern.None,
          lineThickness=1,
          fillColor={28,108,200},
          fillPattern=FillPattern.HorizontalCylinder,
          lineColor={0,0,0}),
        Rectangle(
          extent={{-100,40},{100,30}},
          lineThickness=0.5,
          fillColor={238,46,47},
          fillPattern=FillPattern.Forward,
          lineColor={0,0,0}),
        Rectangle(
          extent={{-100,-30},{100,-40}},
          lineThickness=0.5,
          fillColor={238,46,47},
          fillPattern=FillPattern.Forward,
          lineColor={0,0,0}),
        Rectangle(
          extent={{-100,40},{100,80}},
          lineColor={102,44,145},
          lineThickness=0.5,
          fillColor={0,0,255},
          fillPattern=FillPattern.HorizontalCylinder),
        Rectangle(
          extent={{-100,-80},{100,-40}},
          lineColor={102,44,145},
          lineThickness=0.5,
          fillColor={0,0,255},
          fillPattern=FillPattern.HorizontalCylinder),
        Polygon(
          points={{38,8},{56,0},{38,-8},{38,8}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-44,1},{46,1}},
          color={0,0,0}),
        Line(
          points={{-32,-59},{58,-59}},
          color={28,108,200}),
        Polygon(
          points={{-24,-52},{-42,-60},{-24,-68},{-24,-52}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-32,59},{58,59}},
          color={28,108,200}),
        Polygon(
          points={{-24,66},{-42,58},{-24,50},{-24,66}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid)}));
end ShellTubeHX;
