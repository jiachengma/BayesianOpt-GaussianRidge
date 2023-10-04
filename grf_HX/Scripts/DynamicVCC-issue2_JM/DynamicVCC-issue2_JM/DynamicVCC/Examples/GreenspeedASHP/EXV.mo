within DynamicVCC.Examples.GreenspeedASHP;
model EXV
  import Modelica.Fluid.Utilities.regRoot2;

  extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialIsenthalpicValve;

  Modelica.Blocks.Interfaces.RealInput opening(min=0.0,max=1.0,start=opening_init) "Valve opening" annotation (Placement(transformation(extent={{36,-14},
            {56,6}}),   iconTransformation(extent={{11,-12},{-11,12}},
        rotation=-90,
        origin={1,-46})));
  parameter Real opening_init=0.5;

protected
  Real p_ratio(start=5);

equation
  //Cd=-7.88809953*opening-1.57800815*p_ratio+0.18635886*p_ratio^2+5.59434918;
  Cd=1.8;
  p_ratio=port_a.p/port_b.p;
  m_flow=Av*Cd*opening*regRoot2(2*dp,dp_nominal*1e-9,Medium.density(state_a),Medium.density(state_b));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end EXV;
