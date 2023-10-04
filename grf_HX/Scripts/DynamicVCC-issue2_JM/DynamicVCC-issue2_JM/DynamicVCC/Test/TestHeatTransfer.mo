within DynamicVCC.Test;
model TestHeatTransfer "Test refrigerant heat transfer correlation"
  extends Modelica.Icons.Example;

  replaceable package Medium=DynamicVCC.Media.CoolProp.R410a;

  replaceable model HeatTransfer=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.CorrelationPhaseChange (
  redeclare final model TwoPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Boiling_GungorWinterton);

  Medium.ThermodynamicState states;
  SI.Velocity v;

  HeatTransfer heatTransfer(
  redeclare package Medium=Medium,
  n=1,
  states={states},
  surfaceAreas={0.933},
  dimensions={0.075},
  lengths={4.94},
  vs={v});

  // constant wall temperature heat source
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature wall(T=283);

equation

  states=Medium.setState_ph(10e5,(2+0.25*time)*1e5);

  v=2;

  connect(heatTransfer.heatPorts[1],wall.port);

  annotation (experiment(StopTime=5, __Dymola_Algorithm="Dassl"));
end TestHeatTransfer;
