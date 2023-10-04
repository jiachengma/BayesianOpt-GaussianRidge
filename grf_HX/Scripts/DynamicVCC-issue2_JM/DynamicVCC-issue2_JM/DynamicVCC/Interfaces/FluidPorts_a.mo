within DynamicVCC.Interfaces;
connector FluidPorts_a
  extends DynamicVCC.Interfaces.FluidPort;
  annotation(Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-50,-200},{50,200}},
        initialScale=0.2), graphics={
        Rectangle(
          extent={{50,-200},{-50,200}},
          lineColor={0,127,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-50,180},{50,80}},
          fillColor={0,127,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-50,50},{50,-50}},
          fillColor={0,127,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-50,-80},{50,-180}},
          fillColor={0,127,255},
          fillPattern=FillPattern.Solid)}));
end FluidPorts_a;
