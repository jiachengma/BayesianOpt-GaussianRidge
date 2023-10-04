within DynamicVCC.Test;
model TestPipingModel
  extends Modelica.Icons.Example;

  inner DynamicVCC.Components.System system(
  p_max=45e5,
  p_min=2e5,
  h_max=4.7e5,
  h_min=1.11e5,
  T_max=340,
  T_min=250,
  m_flow_init=0.055,
  m_flow_nominal=0.055,
  momentumDynamics=DynamicVCC.Components.Types.Dynamics.SteadyState_init,
  enableReverseFlow=true);

  replaceable package Medium=DynamicVCC.Media.R410a_NN;

  replaceable model HeatTransfer_1_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    alpha0=500);

  replaceable model HeatTransfer_2_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
   final alpha0=5) "Connect piping heat loss to ambient";

  replaceable model FlowModel_1=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.DetailedFlow;

  Components.Units.HX.Piping pipe(
    redeclare final package Medium_1=Medium,
    final Ncell=5,
    redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
    redeclare final model FlowModel_1=FlowModel_1,
    final modelStructure=DynamicVCC.Components.Types.ModelStructure.av_vb,
    final differentialState=DynamicVCC.Components.Types.DifferentialState.pdh,
    d_o=0.009525,
    d_i=0.0078994,
    length=6.1,
    T_amb=300,
    SteadyState_init=true,
    m_flow_init=0.055,
    p_init=fill(24e5, 5),
    h_init=fill(2.6e5,5),
    Tt_init=fill(310, 5));

    Modelica.Fluid.Sources.MassFlowSource_h source(
    redeclare package Medium=Medium,
    nPorts=1,
    use_m_flow_in=true,
    use_h_in=true);

    Modelica.Fluid.Sources.MassFlowSource_h sink(
    redeclare package Medium=Medium,
    nPorts=1,
    use_m_flow_in=true);
equation
  connect(source.ports[1],pipe.port_a1);
  connect(pipe.port_b1,sink.ports[1]);

  source.m_flow_in=0.055;
  source.h_in=2.6e5;
  sink.m_flow_in=-0.055;




  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StartTime=1,
      StopTime=100,
      Tolerance=0.001,
      __Dymola_Algorithm="Lsodar"));
end TestPipingModel;
