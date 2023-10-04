within DynamicVCC.Interfaces;
connector FluidPort_a
  extends DynamicVCC.Interfaces.FluidPort;
  annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics={Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={0,127,255},
          fillColor={0,127,255},
          fillPattern=FillPattern.Solid), Ellipse(
          extent={{-100,100},{100,-100}},
          fillColor={0,127,255},
          fillPattern=FillPattern.Solid)}));
end FluidPort_a;
