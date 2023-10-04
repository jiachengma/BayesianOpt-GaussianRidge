within DynamicVCC.Test.JCI;
model EXV
  import Modelica.Fluid.Utilities.regRoot2;
  import Modelica.Fluid.Utilities.regPow;
  extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialIsenthalpicValve;

  Modelica.Blocks.Interfaces.RealInput opening(min=0.0,max=1.0) "Valve opening" annotation (Placement(transformation(extent={{36,-14},
            {56,6}}),   iconTransformation(extent={{11,-12},{-11,12}},
        rotation=-90,
        origin={1,-46})));
  parameter SI.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature "Critical temperature";

protected
  SI.TemperatureDifference T_sc_a "subcooling at port_a";

equation
  T_sc_a=Medium.saturationTemperature(Medium.pressure(state_a))-Medium.temperature(state_a);

  /* Discharge coefficient */
  Cd=6.8924*opening-18.3523*opening.^2+15.5418*opening.^3-0.1249;
  //Cd=0.676;
  //Cd=0.421945683435513*opening^0.00172966768*regPow(T_sc_a/Tcrit,-0.1168,1e-4);

  /* Calculate mass flow */
  m_flow=homotopy(
  actual=Av*Cd*opening*regRoot2(2*dp,dp_nominal*1e-5,Medium.density(state_a),Medium.density(state_b)),
  simplified=m_flow_nominal*dp/dp_nominal);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end EXV;
