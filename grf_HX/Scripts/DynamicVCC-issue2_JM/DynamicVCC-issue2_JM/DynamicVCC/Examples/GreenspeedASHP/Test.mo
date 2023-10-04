within DynamicVCC.Examples.GreenspeedASHP;
package Test
  extends Modelica.Icons.ExamplesPackage;
  model TestHX
    extends Modelica.Icons.Example;

    extends DynamicVCC.Examples.GreenspeedASHP.HeatExchangers;

    Modelica.Fluid.Sources.MassFlowSource_h source(
    redeclare package Medium=Medium_1,
    nPorts=1,
    use_m_flow_in=true,
    use_h_in=true);

    Modelica.Fluid.Sources.MassFlowSource_h sink(
    redeclare package Medium=Medium_1,
    nPorts=1,
    use_m_flow_in=true);

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
    connect(source.ports[1], indoorCoil.port_a1);
    connect(indoorCoil.port_b1, sink.ports[1]);
    connect(source_air.ports[1], indoorCoil.ports_a2);
    connect(indoorCoil.ports_b2, sink_air.ports[1]);
    connect(indoorCoil.port_b1, superheat.port);

     source.m_flow_in=0.05;
     //source.h_in=246908;
     source.h_in=442520;
     sink.m_flow_in=-0.05;

    for i in 1:Ncell loop
      //source_air[i].m_flow_in=1.74/Ncell "Outdoor coil";
      //source_air[i].T_in=293.15; // Outdoor coil
      source_air[i].m_flow_in=0.43456/Ncell;
      source_air[i].T_in=295.4;
      source_air[i].X_in={1e-4,1-1e-4};
    end for;

    annotation (experiment(
        StartTime=1,
        StopTime=100,
        __Dymola_Algorithm="Dassl"));
  end TestHX;
end Test;
