within DynamicVCC.Interfaces;
connector HeatPorts_a
  extends HeatPort;
  annotation (Icon(coordinateSystem(extent={{-160,-100},{160,100}}), graphics={
        Rectangle(
          extent={{-205,54},{196,-46}},
          lineColor={127,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-175,49},{-87,-41}},
          lineColor={127,0,0},
          fillColor={127,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-49,49},{39,-41}},
          lineColor={127,0,0},
          fillColor={127,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{78,49},{166,-41}},
          lineColor={127,0,0},
          fillColor={127,0,0},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(extent={{-160,-100},{160,100}})));
end HeatPorts_a;
