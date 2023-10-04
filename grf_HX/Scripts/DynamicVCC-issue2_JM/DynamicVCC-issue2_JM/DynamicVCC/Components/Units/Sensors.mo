within DynamicVCC.Components.Units;
package Sensors
  extends Modelica.Icons.SensorsPackage;
  model T_Superheat
    extends Modelica.Fluid.Sensors.BaseClasses.PartialAbsoluteSensor(redeclare package Medium =
           Modelica.Media.Interfaces.PartialTwoPhaseMedium);

    Modelica.Blocks.Interfaces.RealOutput T(final quantity="ThermodynamicTemperature",
                                            final unit = "K")
      "superheating degree" annotation (Placement(transformation(extent={{60,-10},{80,10}})));

  protected
    Medium.Temperature T_sat;
  equation
    T_sat=Medium.saturationTemperature(port.p);
    T=Medium.temperature(Medium.setState_ph(port.p,inStream(port.h_outflow)))-T_sat;

    annotation (defaultComponentName="temperature",
      Documentation(info="<html>
<p>
This component monitors the temperature of the fluid passing its port.
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</html>"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Line(points={{0,-70},{0,-100}}, color={0,0,127}),
          Ellipse(
            extent={{-20,-98},{20,-60}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-12,40},{12,-68}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,40},{-12,80},{-10,86},{-6,88},{0,90},{6,88},{10,86},{
                12,80},{12,40},{-12,40}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Line(
            points={{-12,40},{-12,-64}},
            thickness=0.5),
          Line(
            points={{12,40},{12,-64}},
            thickness=0.5),
          Line(points={{-40,-20},{-12,-20}}),
          Line(points={{-40,20},{-12,20}}),
          Line(points={{-40,60},{-12,60}}),
          Line(points={{12,0},{60,0}}, color={0,0,127})}),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Ellipse(
            extent={{-20,-88},{20,-50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-12,50},{12,-58}},
            lineColor={191,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,50},{-12,90},{-10,96},{-6,98},{0,100},{6,98},{10,96},{
                12,90},{12,50},{-12,50}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Line(
            points={{-12,50},{-12,-54}},
            thickness=0.5),
          Line(
            points={{12,50},{12,-54}},
            thickness=0.5),
          Line(points={{-40,-10},{-12,-10}}),
          Line(points={{-40,30},{-12,30}}),
          Line(points={{-40,70},{-12,70}}),
          Text(
            extent={{102,-16},{-18,-46}},
            lineColor={0,0,0},
            textString="SH"),
          Text(
            extent={{-150,106},{150,146}},
            textString="%name",
            lineColor={0,0,255}),
          Line(points={{12,0},{60,0}}, color={0,0,127})}),
                Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end T_Superheat;

  model Temperature_grid "Average temperature of thermocouple grid"
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
      "Medium in the sensor";

    parameter Integer nPorts=1;

    DynamicVCC.Interfaces.FluidPorts_a ports[nPorts](redeclare each package Medium = Medium) annotation (Placement(transformation(
          origin={-8,-100},
          extent={{-10,-10},{10,10}},
          rotation=90), iconTransformation(
          extent={{-6,-18},{6,18}},
          rotation=90,
          origin={0,-96})));

    Modelica.Blocks.Interfaces.RealOutput T(final quantity="ThermodynamicTemperature",
                                            final unit = "K", min=0)
                                            annotation (Placement(transformation(extent={{60,-10},{80,10}})));

  protected
    Medium.Temperature Tport[nPorts];

  equation
    for i in 1:nPorts loop
      ports[i].m_flow = 0;
      ports[i].h_outflow = Medium.h_default;
      ports[i].Xi_outflow = Medium.X_default[1:Medium.nXi];
      ports[i].C_outflow = zeros(Medium.nC);
      Tport[i] = Medium.temperature(Medium.setState_phX(ports[i].p, inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)));
    end for;

    T=sum(Tport)/nPorts;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Ellipse(
            extent={{-20,-88},{20,-50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={217,67,180},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-12,50},{12,-58}},
            lineColor={191,0,0},
            fillColor={217,67,180},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-12,50},{-12,90},{-10,96},{-6,98},{0,100},{6,98},{10,96},{
                12,90},{12,50},{-12,50}},
            lineColor={0,0,0},
            lineThickness=0.5),
          Line(
            points={{-12,50},{-12,-54}},
            thickness=0.5),
          Line(
            points={{12,50},{12,-54}},
            thickness=0.5),
          Line(points={{-40,-10},{-12,-10}}),
          Line(points={{-40,30},{-12,30}}),
          Line(points={{-40,70},{-12,70}}),
          Text(
            extent={{102,-16},{-18,-46}},
            lineColor={0,0,0},
            textString="T"),
          Line(points={{12,0},{60,0}}, color={0,0,127}),
          Text(
            extent={{-150,102},{150,142}},
            textString="%name",
            lineColor={0,0,255})}),                                Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Temperature_grid;
end Sensors;
