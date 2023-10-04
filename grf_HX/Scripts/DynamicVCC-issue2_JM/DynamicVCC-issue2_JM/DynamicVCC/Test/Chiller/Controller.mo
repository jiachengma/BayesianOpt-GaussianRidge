within DynamicVCC.Test.Chiller;
model Controller
 extends Modelica.Blocks.Interfaces.SI2SO; //u1=Twater, u2=Tset, y=gamma

 parameter Real starttime=1;
 parameter Real gamma_init;

protected
  parameter Real gamma_max=1;
  parameter Real gamma_min=0.05;
  Real com_action;
  Real act_pre;
  Real gamma(start=gamma_init);
  Real error;

  Real Tw;
  Real Tset;
  Boolean open(start=true);
  parameter Integer sampletime_act=15;
  parameter Integer sampletime_vane=1;
algorithm
  Tw:=u1;
  Tset:=u2;
  when sample(starttime,sampletime_act) then
    error:=Tw - Tset;
    open:=if error > 0 then true else false;
    act_pre:=vaneAction(abs(error));
  end when;

  when sample(starttime,sampletime_vane) then
    (gamma,com_action):=controlSignal(open,gamma,com_action);
  end when;

  y:=max(min(gamma_max,gamma),gamma_min);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Controller;
