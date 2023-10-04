within DynamicVCC.Interfaces;
connector HeatPorts_b
  extends HeatPort;
  annotation (Icon(coordinateSystem(extent={{-160,-100},{160,100}}), graphics={
        Rectangle(
          extent={{-198,50},{202,-51}},
          lineColor={127,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-168,44},{-80,-46}},
          lineColor={127,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-42,46},{46,-44}},
          lineColor={127,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{84,45},{172,-45}},
          lineColor={127,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(extent={{-160,-100},{160,100}})));
end HeatPorts_b;
