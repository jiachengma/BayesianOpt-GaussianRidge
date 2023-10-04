within DynamicVCC.Interfaces;
connector HeatPort_in
  extends DynamicVCC.Interfaces.HeatPort;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}), graphics={Rectangle(
          extent={{-50,50},{50,-50}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid), Text(
          extent={{-120,120},{100,60}},
          textColor={191,0,0},
          textString="%name")}));
end HeatPort_in;
