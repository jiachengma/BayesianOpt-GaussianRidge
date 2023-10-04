within DynamicVCC.Test.Chiller;
model Bulb "R500 bulb"
  extends Modelica.Blocks.Interfaces.SISO;

  parameter Real C=100;
  parameter SI.Temperature T_init=280;
  parameter Boolean SteadyState_init=false;
  SI.Temperature T(start=T_init);
equation
  der(T) = (u - T) / C;
  y =( 0.00085 * T^3 - 0.5466*T^2 + 120.3*T - 9038)*1e3; //R500

initial equation
  if SteadyState_init then
    der(T)=0;
  else
    T=T_init;
  end if;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Bulb;
