within DynamicVCC.Test.JCI.Test.TPWL;
model LinearizeHX "Linearize heat exchanger and pipe model"
  extends Modelica.Icons.Example;

  //Modelica.Blocks.Interfaces.RealInput u[5](start={0.05,22.8916e5,4.75e5,2,305});
  Modelica.Blocks.Interfaces.RealInput u[5](start={0,19e5,2.4e5,0,305});

  //Modelica.Blocks.Interfaces.RealInput u[5](start={0.02364,2.31e6,4.75e5,2,305});
  Modelica.Blocks.Interfaces.RealOutput y[3];

  extends DynamicVCC.Test.JCI.Test.TPWL.HeatExchangers(
    final p_a_start_OD=p_a_start_OD_ss,
    final p_b_start_OD=p_b_start_OD_ss,
    final h_init_OD=h_init_OD_ss,
    final Tt_init_OD=Tt_init_OD_ss,
    final p_a_start_ID=p_a_start_ID_ss,
    final p_b_start_ID=p_b_start_ID_ss,
    final h_init_ID=h_init_ID_ss,
    final Tt_init_ID=Tt_init_ID_ss);

  // Initial conditions
/*
  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=2.3576e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=2.30503e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={447244.875000000,432615.343750000,420101.875000000,407676.781250000,395340.218750000,383092.187500000,370932.843750000,358862.250000000,346880.468750000,334987.656250000,323183.812500000,311469.093750000,299843.531250000,288307.250000000,276860.343750000,265502.875000000,259624.609375000,256170.562500000,254155.765625000,252985.593750000};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={319.686248779297,312.780456542969,311.655059814453,311.608032226563,311.560974121094,311.513854980469,311.466705322266,311.419525146484,311.372283935547,311.324981689453,311.277648925781,311.230255126953,311.182830810547,311.135345458984,311.087829589844,311.040252685547,308.126251220703,306.836975097656,306.071533203125,305.622344970703};
  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=9.6703e5;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=9.4187e5;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={266489.250000000,277491.781250000,288509.781250000,299545.500000000,310601.125000000,321678.937500000,332781.218750000,343910.187500000,355068.218750000,366257.562500000,377480.625000000,388739.718750000,400037.281250000,411375.718750000,422757.468750000,427689.718750000,431256.906250000,433832.937500000,435692.843750000,437037.375000000};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={280.202697753906,280.185455322266,280.165313720703,280.142272949219,280.116333007813,280.087493896484,280.055694580078,280.020904541016,279.983154296875,279.942382812500,279.898529052734,279.851623535156,279.801605224609,279.748413085938,279.692047119141,288.083007812500,289.859039306641,291.148529052734,292.080200195313,292.750701904297};
  */

  parameter Medium_1.AbsolutePressure p_a_start_OD_ss=1.9833e6;
  parameter Medium_1.AbsolutePressure p_b_start_OD_ss=1.9834e6;
  parameter Medium_1.SpecificEnthalpy h_init_OD_ss[Ncell]={426388.250000000,426388.250000000,426388.250000000,426388.218750000,426388.218750000,426388.218750000,426388.218750000,426388.187500000,426388.187500000,426388.187500000,426388.156250000,426388.156250000,426388.156250000,426388.125000000,426388.125000000,426388.093750000,426388.093750000,426388.062500000,426388.031250000,337495.593750000};
  parameter SI.Temperature Tt_init_OD_ss[Ncell]={305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,304.999969482422,304.999969482422};

  parameter Medium_1.AbsolutePressure p_a_start_ID_ss=1.9834e6;
  parameter Medium_1.AbsolutePressure p_b_start_ID_ss=1.9834e6;
  parameter Medium_1.SpecificEnthalpy h_init_ID_ss[Ncell]={255435.343750000,254673.437500000,254366.328125000,254086.750000000,253857.171875000,253722.187500000,253546.343750000,253381.296875000,253226.093750000,253077.140625000,252939.312500000,252811.390625000,252693.312500000,252589.906250000,252504,252442.500000000,252409.468750000,252435.328125000,252488.515625000,252368.078125000};
  parameter SI.Temperature Tt_init_ID_ss[Ncell]={305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305,305};
  parameter Medium_1.MassFlowRate m_flow_init_ss=0;

  // Connect piping between indoor and outdoor units
  replaceable model HeatTransfer_1_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    alpha0=50);

  replaceable model HeatTransfer_2_piping=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
   final alpha0=5) "Connect piping heat loss to ambient";

  replaceable model FlowModel_1_LL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
  final dp_nominal=1.3796e4/Ncell_piping,
  final k=1,
  final b=2);

  replaceable model FlowModel_1_VL=DynamicVCC.Components.Pipes.BaseClasses.FlowModels.NominalFrictionFlow (
  final dp_nominal=2.0953e4/Ncell_piping,
  final k=1,
  final b=2);

  parameter Integer Ncell_piping=3;

  Components.Units.HX.Piping pipe(
    redeclare final package Medium_1=Medium_1,
    final Ncell=Ncell_piping,
    redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
    redeclare final model FlowModel_1=FlowModel_1_LL,
    redeclare final model SlipRatio=SlipRatio,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final useLumpedPressure=useLumpedPressure,
    d_o=0.009525,
    d_i=0.0078994,
    length=6.1,
    T_amb=300,
    p_a_start=p_b_start_OD,
    h_init=fill(h_init_OD[Ncell],Ncell_piping),
    Tt_init=fill(Tt_init_OD[Ncell], Ncell_piping));
/*
  Components.Units.HX.Piping pipe(
  redeclare final package Medium_1=Medium_1,
    final Ncell=Ncell_piping,
    redeclare final model HeatTransfer_1=HeatTransfer_1_piping,
    redeclare final model HeatTransfer_2=HeatTransfer_2_piping,
    redeclare final model FlowModel_1=FlowModel_1_VL,
    redeclare final model SlipRatio=SlipRatio,
    final modelStructure=modelStructure,
    final differentialState=differentialState,
    final useLumpedPressure=useLumpedPressure,
    d_o=0.009525,
    d_i=0.0078994,
    length=6.1,
    T_amb=300,
    p_a_start=p_b_start_ID,
    h_init=fill(h_init_ID[Ncell],Ncell_piping),
    Tt_init=fill(Tt_init_ID[Ncell], Ncell_piping));
*/

  Modelica.Fluid.Sources.MassFlowSource_h source(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_m_flow_in=true,
  use_h_in=true);
/*
  Modelica.Fluid.Sources.MassFlowSource_h sink(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_m_flow_in=true);
*/

  Modelica.Fluid.Sources.Boundary_ph sink(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_p_in=true);

  Modelica.Fluid.Sources.MassFlowSource_T source_air[Ncell](
  redeclare each final package Medium=Medium_2,
  each use_m_flow_in=true,
  each use_T_in=true,
  each use_X_in=true,
  each nPorts=1);

  Modelica.Fluid.Sources.Boundary_pT sink_air[Ncell](
  redeclare each package Medium=Medium_2,
  each nPorts=1);

  DynamicVCC.Components.Units.Sensors.T_Superheat superheat(
  redeclare package Medium=Medium_1);

equation

  connect(source.ports[1],HX.port_a1);
  connect(HX.port_b1,pipe.port_a1);
  connect(pipe.port_b1,sink.ports[1]);
  connect(source_air.ports[1],HX.ports_a2);
  connect(HX.ports_b2,sink_air.ports[1]);
  connect(HX.port_b1,superheat.port);

/*
  source.m_flow_in=0.05;
  source.h_in=4.75e5;
  //source.h_in=2.555e5;
  sink.m_flow_in=-0.05;
  //sink.h_in=2.55e5;
  //sink.p_in=23.3e5;
  //sink.p_in=9.4e5;
  */

 source.m_flow_in=u[1];
 //sink.m_flow_in=u[2];
 sink.p_in=u[2];
 source.h_in=u[3];

  for i in 1:Ncell loop
    source_air[i].m_flow_in=u[4]/Ncell "Outdoor coil";
    //source_air[i].m_flow_in=1/Ncell "Indoor coil";
    source_air[i].T_in=u[5]; // Outdoor coil
    //source_air[i].T_in=294.5; // Indoor coil
    source_air[i].X_in={1e-4,1-1e-4};
  end for;

  y[1]=HX.port_a1.p;
  y[2]=HX.Ta_out_ave;
  y[3]=superheat.T;

  annotation (experiment(
      StartTime=1,
      StopTime=3000,
      __Dymola_Algorithm="Dassl"));
end LinearizeHX;
