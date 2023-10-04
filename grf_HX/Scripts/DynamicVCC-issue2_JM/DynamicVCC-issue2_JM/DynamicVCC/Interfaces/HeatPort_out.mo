within DynamicVCC.Interfaces;
connector HeatPort_out
  extends DynamicVCC.Interfaces.HeatPort;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={191,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end HeatPort_out;
