within DynamicVCC.Components.Units;
package SteadyMassFlowPressure "Static component models using steady mass flow pressure"
  partial model SISOFlow "Base class for one-dimensional component"

    extends DynamicVCC.Interfaces.PartialTwoPort;

    parameter Real L=1 "Inertance";

    parameter Medium.MassFlowRate m_flow_init=1;

    parameter Boolean SteadyState_init=false;

  protected
    SI.Pressure r;
    Medium.MassFlowRate m_flow(min=0);
    Medium.AbsolutePressure p_in "Steady mass flow pressure";
    Medium.AbsolutePressure p_out "Steady mass flow pressure";
    Medium.SpecificEnthalpy h_in;
    Medium.SpecificEnthalpy h_out;
    SI.Pressure dp;

  equation

    port_a.m_flow+port_b.m_flow=0;

    r + L*der(m_flow) = 0;

    p_out = p_in + dp;

    r = port_b.p - p_out;

    // Boundary conditions
    p_in=port_a.p;
    h_in=inStream(port_a.h_outflow);
    port_a.m_flow=m_flow;

  initial equation
    if SteadyState_init then
      der(m_flow)=0;
    else
      m_flow=m_flow_init;
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end SISOFlow;

  model CheckValve
    extends SISOFlow;

    Modelica.Blocks.Interfaces.RealInput u(max=1.0,min=0.0);

    parameter SI.Area A_open=1 "Fully opened area";
    parameter Real opening_min=1e-4;

  protected
    Medium.ThermodynamicState state_a;
    Medium.Density rho=Medium.density(state_a);
    Real opening(max=1.0,min=opening_min);

  equation
    u=opening_min+(1.0-opening_min)*opening;

    state_a=Medium.setState_ph(p_in,h_in);

    dp=-(m_flow/(u*A_open))^2/rho;

    h_out=h_in;

    // Boundary conditions
    port_b.h_outflow=h_out;
    port_a.h_outflow=Medium.h_default;

  end CheckValve;

  model TestValve
    extends Modelica.Icons.Example;

    package Medium=DynamicVCC.Media.CoolProp.R410a;

    Modelica.Fluid.Sources.Boundary_ph source(
    nPorts=1,
    redeclare package Medium=Medium,
    use_p_in=true,
    use_h_in=true);

    Modelica.Fluid.Sources.Boundary_ph sink(
    nPorts=1,
    redeclare package Medium=Medium,
    use_p_in=true);

    parameter SI.Area Av=4e-6;

    DynamicVCC.Components.Units.MassFlowDevices.Valve.CheckValve_Rev valve(redeclare package Medium = Medium, A_open=Av);

  equation
    valve.u=1.0;
    source.p_in=32e5;
    source.h_in=4.2e5;
    sink.p_in=25e5;

    connect(source.ports[1],valve.port_a);
    connect(sink.ports[1],valve.port_b);

    annotation (experiment(StopTime=10, __Dymola_Algorithm="Dassl"));
  end TestValve;

  model ReversingValve

    extends DynamicVCC.Interfaces.PartialTwoPort;

    Modelica.Blocks.Interfaces.RealInput opening_SP(max=1,min=0) annotation (Placement(transformation(extent={{-82,56},{-62,76}}),  iconTransformation(extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-48,86})));

    TransientVCC.Interfaces.FluidPort_a port_a2(redeclare package Medium=Medium)
    "Inlet of the second pair of ports, connect to outdoor coil" annotation (Placement(transformation(extent={{90,52},{110,72}}),  iconTransformation(extent={{90,52},{110,72}})));

    TransientVCC.Interfaces.FluidPort_b port_b2(redeclare package Medium=Medium)
    "Outlet of the first pair of ports, connect to indoor coil" annotation (Placement(transformation(extent={{90,-78},{110,-58}}),  iconTransformation(extent={{90,-78},{110,-58}})));

    parameter SI.Area A_open "Area of each valve when fully opened";
    parameter SI.Time tau "Time constant of first-order filter";
    parameter Real L "Inertance";
    parameter Real opening_min=1e-4;

    parameter Medium.MassFlowRate m_flow_init=1;
    parameter Real opening_init=1;
    parameter Boolean SteadyState_init=false;

  protected
    Real opening(max=1.0,min=0.0,start=opening_init);

  public
    DynamicVCC.Components.Units.SteadyMassFlowPressure.CheckValve checkValve_dis[2](
      redeclare each package Medium = Medium,
      each SteadyState_init=SteadyState_init,
      each A_open=A_open,
      each L=L,
      each opening_min=opening_min,
      m_flow_init={m_flow_init,1e-5});

    DynamicVCC.Components.Units.SteadyMassFlowPressure.CheckValve checkValve_suc[2](
      redeclare each package Medium = Medium,
      each SteadyState_init=SteadyState_init,
      each A_open=A_open,
      each L=L,
      each opening_min=opening_min,
      m_flow_init={1e-5,m_flow_init});

  equation

    der(opening)=(opening_SP-opening)/tau;

    checkValve_dis[1].opening=opening;
    checkValve_dis[2].opening=1.0-opening;
    checkValve_suc[1].opening=1.0-opening;
    checkValve_suc[2].opening=opening;

    // Boundary conditions
    connect(port_a,checkValve_dis[1].port_a);
    connect(port_a,checkValve_dis[2].port_a);
    connect(port_b,checkValve_suc[1].port_b);
    connect(port_b,checkValve_suc[2].port_b);
    connect(port_a2,checkValve_dis[2].port_b);
    connect(port_a2,checkValve_suc[2].port_a);
    connect(port_b2,checkValve_dis[1].port_b);
    connect(port_b2,checkValve_suc[1].port_a);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-90,-22},{-10,18}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-10,78},{30,-82}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{30,78},{90,50}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{30,12},{90,-16}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{30,-54},{90,-82}},
            lineColor={102,44,145},
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid)}),                      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ReversingValve;

  model BasicCheckValve

    extends DynamicVCC.Components.Units.SteadyMassFlowPressure.SISOFlow;

    parameter SI.Area A_open=1 "Fully opened area";
    parameter Real opening_min=1e-4;

    parameter Medium.MassFlowRate m_flow_ref=0.037;
    parameter SI.Pressure dp_ref=0.2e5;
    parameter Medium.Density rho_ref=100;

  protected
    Medium.ThermodynamicState state_a;
    Medium.Density rho=Medium.density(state_a);

  equation
    state_a=Medium.setState_ph(p_in,h_in);

    dp=min(0.0,-(m_flow/((opening_min+(1.0-opening_min)*opening)*A_open))^2/rho);

    h_out=h_in;

    // Boundary conditions
    port_b.h_outflow=h_out;
    port_a.h_outflow=Medium.h_default;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end BasicCheckValve;

  model Compressor

    extends SISOFlow;

    import Modelica.Math.Vectors;

    // Compressor speed
    Modelica.Blocks.Interfaces.RealInput speed(quantity="Frequency",final unit="Hz",start=speed_init) annotation (Placement(transformation(extent={{54,14},
            {74,34}}),
            iconTransformation(extent={{90,-86},{76,-72}})));

    input SI.Temperature T_amb;

    parameter SI.Volume Vs "Displacement";

    // Initialization
    parameter SI.Power Pwr_init=1000;
    parameter SI.Frequency speed_init=53;

  /******************* Maps *********************/
    constant Real MASS_map[:,:]=[159.291276008480,1.78973432192720,-2.63300401447600,-0.00783168526800000,0.00655027705640000,0.0237112803532000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05;202.524738043980,2.38662488486400,-2.52904562081600,0.000788056896000000,0.00982541558460000,0.0219667500162000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05;208.227524860960,3.90368191664000,-1.88446018172060,0.0166412097874000,-0.0105906570088000,0.0193042983958000,0.000299828196800000,-7.68342117400000e-05,0.000106149940600000,-8.25809736400000e-05;276.874273089700,3.58040601073760,-2.19637876110400,0.0180275432788000,0.0163756926410000,0.0184776903696000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05;353.334649209520,5.23843535200040,-1.45840248837800,0.0419712739094000,0.0254732990930000,0.0136317733376000,0.000304540880600000,0,-3.04829580000000e-05,-7.37570460000000e-05];
    constant Real POWER_map[:,:]=[-279.683678251600,2.10990406957400,13.8321489984790,-0.109808251721000,0.0201255194287000,-0.0794447884860000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;-24.4061357146600,-0.108011374730000,12.4036190197150,-0.109808251721000,0.0301882786604000,-0.0524247383310000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;198.535720623520,-1.76434550257700,11.2779341447780,-0.109808251721000,0.0395242837934000,-0.0273561290030000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;448.101327449020,-3.19728107192900,10.2163985612780,-0.109808251721000,0.0503137980891000,0.00161536680550000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000;1122.45682679060,-4.50942876811700,8.66023324954760,-0.109808251721000,0.0782659091039000,0.0766710684376000,0,0.000596006144300000,-0.000461548542000000,0.000588949801300000];
    constant SI.Frequency frequency[:]={0,30,45,60,75,117};
    parameter Real Fcor=0.75;

  protected
    Medium.MassFlowRate m_dot_map[size(frequency,1)];
    SI.Power Pwr_map[size(frequency,1)];
    SI.Power Pwr;
    SI.HeatFlowRate Q_loss "Heat loss";
    Medium.SaturationProperties sat_suc=Medium.setSat_p(p_in);
    Medium.ThermodynamicState state_suc=Medium.setState_ph(p_in,h_in);
    Medium.SaturationProperties sat_dis=Medium.setSat_p(p_out);
    //Medium.ThermodynamicState state_dis=Medium.setState_ph(p_out,h_out);

    Medium.Temperature Te_sat=Medium.saturationTemperature_sat(sat_suc);
    Medium.Temperature Tc_sat=Medium.saturationTemperature_sat(sat_dis);
    Medium.Density rho_suc=Medium.density(state_suc);
    Real p_ratio(min=0) "Ratio of discharge and suction pressures";
    Real f_Qloss(start=0.55) "Heat loss coefficient";
    //Real UA "Chamber outter wall heat transfer coefficient";
    Real eta_v(start=0.8);
  equation

    m_dot_map = TransientVCC.Utilities.map10coef(
        MASS_map,
        Te_sat,
        Tc_sat) ./ 7937;
    Pwr_map = TransientVCC.Utilities.map10coef(
        POWER_map,
        Te_sat,
        Tc_sat);

    eta_v=-0.0308*p_ratio+0.0041*p_ratio^2-0.0873*speed/53+0.9553;
    eta_v=min(1.0,max(0.5,m_flow/(Vs*Medium.density(state_suc).*speed)));

    Pwr=Vectors.interpolate(frequency,Pwr_map,speed);

    p_out=p_ratio*p_in;
    //f_Qloss=0.2587*m_flow/0.035-0.4677*T_amb/273+0.6992;
    //f_Qloss=min(1.0,max(0.0,-0.6214*m_flow/0.035+9.2704*T_amb/344.494+0.6276*speed/53-7.0899));
    f_Qloss=0.56;
    Q_loss=Pwr*f_Qloss;
    Pwr=Q_loss+m_flow*(h_out-h_in);

    // Boundary conditions
    port_b.h_outflow=h_out;
    port_a.h_outflow=Medium.h_default;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StartTime=10050,
        StopTime=29400,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end Compressor;

  model ExpansionValve

    extends SISOFlow;

    Modelica.Blocks.Interfaces.RealInput u(max=1.0,min=0.0);

    parameter SI.Area A_open=1 "Fully opened area";
    parameter Real opening_min=1e-4;

  protected
    Medium.ThermodynamicState state_a;
    Medium.ThermodynamicState state_b;
    Real opening(max=1.0,min=opening_min);

  equation
    u=opening_min+(1.0-opening_min)*opening;

    state_a=Medium.setState_ph(p_in,h_in);

    dp=-(m_flow/(u*A_open))^2/rho;

    h_out=h_in;

    // Boundary conditions
    port_b.h_outflow=h_out;
    port_a.h_outflow=Medium.h_default;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ExpansionValve;

  model Fan

    import Modelica.Fluid.Utilities.regStep;
    import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
    replaceable package Medium=Modelica.Media.Air.MoistAir;

    Modelica.Blocks.Interfaces.RealInput T_in(unit="K",min=250,start=271);

    Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](each unit="1");

    Modelica.Blocks.Interfaces.RealInput speed "RPM";

    parameter Integer Ncell=1;
    parameter Medium.AbsolutePressure p_atm=Medium.p_default;
    parameter Medium.Density rho_nominal=1.23;
    parameter Real eta=0.0756 "Fan efficiency";
    parameter Real L=1 "Inertance";

    // Initial conditions
    parameter SI.VolumeFlowRate V_dot_init=1;

    // Connector
    TransientVCC.Interfaces.FluidPorts_b ports[Ncell](redeclare each final package Medium=Medium,
    each m_flow(max=0,start=-V_dot_init*rho_nominal/Ncell));

  protected
    Medium.Density rho "Upstream flow density";
    Medium.SpecificEnthalpy h "Upstream flow enthalpy";
    SI.Pressure p_drop(min=0,start=7);
    SI.VolumeFlowRate V_dot(start=V_dot_init);
    SI.VolumeFlowRate V_dot_nominal(start=V_dot_init);
    SI.VolumeFlowRate V_dot_curve(start=V_dot_init);
    Medium.MassFlowRate m_flows[Ncell](each min=0);
    SI.Pressure r[Ncell],r_mix,r_C,p_mix,p_C(start=p_atm);
    Medium.MassFlowRate m_f(start=V_dot_init*rho_nominal,min=0);
    SI.Power Pwr;

  equation
  /*
  p_drop=ports[1].p-p_atm;
  for i in 1:Ncell-1 loop
    ports[i+1].p=ports[1].p;
  end for;
  */
    for i in 1:Ncell loop
      r[i]+L*der(m_flows[i])=0;
      r[i]+ports[i].p=r_mix+p_mix;
    end for;
    r_C+p_C=p_atm;
    p_drop=p_mix-p_C;
    -L*der(m_f)=r_C-r_mix;
    sum(der(m_flows))=der(m_f);
    m_f=V_dot*rho;
    p_mix=sum(m_flows.*ports.p)/sum(m_flows);

    rho=Medium.density(Medium.setState_pTX(p_atm,T_in,X_in));
    h=Medium.specificEnthalpy(Medium.setState_pTX(p_atm,T_in,X_in));

    Pwr=V_dot*p_drop/eta;
    V_dot=V_dot_nominal;

    V_dot_nominal=0.0027*speed;

    // Fan curve
     V_dot_curve=min(2,max(1e-7,-0.0009489*p_drop^4 + 0.03018*p_drop^3 -0.3424*p_drop^2 + 1.558*p_drop -1.09));

    // Boundary conditions
    for i in 1:Ncell loop
      ports[i].h_outflow=h;
      ports[i].Xi_outflow[Medium.Water]=X_in[Medium.Water];
      ports[i].m_flow=-max(1e-5,m_flows[i]);
    end for;

  initial equation
    m_flows=fill(V_dot_init*rho_nominal/Ncell,Ncell);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(
            preserveAspectRatio=false)),
      experiment(
        StartTime=622,
        StopTime=10000,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end Fan;

  model Fan_new

    import Modelica.Fluid.Utilities.regStep;
    import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
    replaceable package Medium=Modelica.Media.Air.MoistAir;

    Modelica.Blocks.Interfaces.RealInput speed "RPM";

    parameter Integer Ncell=1;
    parameter Medium.AbsolutePressure p_atm=Medium.p_default;
    parameter Medium.Density rho_nominal=1.23;
    parameter Medium.AbsolutePressure p_rise_nominal=6;
    parameter Real eta=0.0756 "Fan efficiency";
    parameter Real L=1 "Inertance";
    parameter Real speed_nominal=492.5 "RPM";

    // Initial conditions
    parameter SI.VolumeFlowRate V_dot_init=1;

    // Connector
    ThermofluidStream.Interfaces.Inlet inlet(redeclare final package Medium=Medium,
    m_flow(min=0,start=V_dot_init*rho_nominal));

    ThermofluidStream.Interfaces.Outlet outlet(redeclare final package Medium=Medium,
    m_flow(max=0,start=-V_dot_init*rho_nominal));

  protected
    Medium.Density rho "Upstream flow density";
    SI.Pressure p_rise(min=0,start=p_rise_nominal);
    SI.VolumeFlowRate V_dot(start=V_dot_init);
    SI.VolumeFlowRate V_dot_nominal(start=V_dot_init);
    SI.VolumeFlowRate V_dot_curve(start=V_dot_init);
    SI.Power Pwr;
    Medium.AbsolutePressure p_out;

  equation

    //Mass balance
    inlet.m_flow+outlet.m_flow=0;

    -L*der(inlet.m_flow)=outlet.r-inlet.r;

    inlet.m_flow=rho*V_dot;

    p_out-Medium.pressure(inlet.state)=regStep(speed-491,p_rise,p_rise_nominal*speed/speed_nominal,1e-3);

    rho=Medium.density(inlet.state);
    Pwr=V_dot*(p_out-Medium.pressure(inlet.state))/eta;
    //V_dot=regStep(speed-491,V_dot_curve,V_dot_nominal,1e-3);
    V_dot_curve=V_dot;
    V_dot_nominal=0.0027*speed;

    // Fan curve
    //V_dot_curve=min(2,max(1e-7,-0.0009489*p_rise^4 + 0.03018*p_rise^3 -0.3424*p_rise^2 + 1.558*p_rise -1.09));
    p_rise = -400*V_dot_curve^3 + 1338*V_dot_curve^2 -1494*V_dot_curve + 564.3;
    // Boundary conditions
    outlet.state=Medium.setState_phX(p_out,Medium.specificEnthalpy(inlet.state),{inlet.state.X[Medium.Water]});

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Fan_new;

  model MoistAirCrossFlow "Moist air heat and mass transfer for cross flow HX"

    import Modelica.Constants.T_zero;
    import Modelica.Fluid.Utilities.regStep;

    parameter Integer Ncell=1 "Number of control volumes";

    // extending dry coil model
    extends DynamicVCC.Components.Units.SteadyMassFlowPressure.BaseClasses.DryAirCoil(final n=Ncell);

    output SI.MassFlowRate m_flows_dehumid[n];

    parameter SI.LewisNumber Le=system.Le "Lewis number";
    parameter SI.SpecificEnthalpy delta_h_ig=2836.6e+3 "Latent heat of sublimation";

    SI.SpecificEnthalpy delta_h_fg[n]=Medium.enthalpyOfVaporization(T_a) "Latent heat of condensation";
    Medium.MassFraction w_s[n];

  protected
      Real Ntu_mass[n];
      SI.HeatFlowRate Q_flows_lat[n] "Latent heat transfer rate";

  equation

  /******************** Mass transfer of moist air ************************/
      for i in 1:n loop
        w_s[i]=Medium.xsaturation_pT(Medium.pressure(states_a[i]), min(T_max,max(T_min,T_s[n])));
      end for;
      m_flows_dehumid={m_flows[i]*(w_a[i]-w_b[i]) for i in 1:n};
      Ntu_mass=Ntu/Le^(2/3);
      w_b={w_a[i]+min(0,(w_s[i]-w_a[i]))*(1-exp(-Ntu_mass[i])) for i in 1:n};

  /************************ Heat transfer **************************/
      Q_flows_lat={m_flows[i]*(w_b[i]-w_a[i])*regStep(T_s[i]+T_zero,delta_h_fg[i],delta_h_ig,1e-2) for i in 1:n};
      Q_flows_forced=Q_flows_sen+Q_flows_lat;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StartTime=1700,
        StopTime=5000,
        __Dymola_Algorithm="Radau"));
  end MoistAirCrossFlow;

  model FinTubeHX "1-D cross-flow fin-and-tube heat exchanger"

    // extending refrigerant flow models and tube wall models
    extends DynamicVCC.Components.Units.HX.BaseClasses.PartialHX(final C_metalWall=C_FinTube);

    // Tube and fin geometry
    parameter SI.Mass M_metalWall "Tube wall mass";
    parameter SI.SpecificHeatCapacity cp_metalWall "Tube material specific heat";
    parameter SI.Area Ac_2 "Air side cross-sectional area";
    parameter Real Eta_fin_overall=1 "Overall fin efficiency";
    parameter SI.Area As_2 "Air side heat transfer area";
    parameter SI.Diameter diameter_2;
    parameter SI.Length L_2;
    parameter SI.Mass M_fin "Fin mass";
    parameter SI.SpecificHeatCapacity cp_fin=1 "Fin material specific heat";
    parameter SI.LewisNumber Le=system.Le "Lewis number";

    //input SI.Area Ac_eff[Ncell];
  protected
    parameter SI.HeatCapacity C_FinTube=M_metalWall*cp_metalWall+M_fin*cp_fin "Heat capacity of tube and fin assuming uniform temperature";

    //Initial conditions
  public
    parameter SI.MassFlowRate m_flows_a_init[Ncell]=ones(Ncell) "Initial air flow rate";
    parameter SI.ThermodynamicTemperature Ta_init[Ncell]=fill(Medium_2.reference_T,Ncell) "Initial air exit temperature";
    parameter SI.MassFraction w_a_init[Ncell]=fill(1e-3,Ncell) "Initial air exit humidity ratio";
    SI.ThermodynamicTemperature Ta_out_ave(start=Ta_init[1]);

    /*  Air side  */

    replaceable model HeatTransfer_2 =
        TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_ConstCoefficient
    "Air side heat transfer model";

    replaceable model Friction_2 =
        TransientVCC.Component.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor
    "Air side frictional pressure drop model";

    // Air flow
    replaceable model AirFlow =
        DynamicVCC.Components.Units.SteadyMassFlowPressure.MoistAirCrossFlow;

    AirFlow airFlow(redeclare final package Medium=Medium_2,
    final Ncell=Ncell,
    redeclare final model HeatTransfer=HeatTransfer_2,
    redeclare final model FrictionalDP=Friction_2,
    final Eta_fin_overall=Eta_fin_overall,
    final Le=Le,
    final diameters=fill(diameter_2,Ncell),
    final length=fill(L_2,Ncell),
    final As=fill(As_2/Ncell,Ncell),
    final Ac=fill(Ac_2/Ncell,Ncell),
    final m_flows_init=m_flows_a_init,
    final T_out_init=Ta_init,
    final w_out_init=w_a_init);

    // Air flow connectors
    ThermofluidStream.Interfaces.Inlet inlets[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
          extent={{-8,-32},{8,32}},
          rotation=90,
          origin={0,-68})));

    ThermofluidStream.Interfaces.Outlet outlets[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
          extent={{-8,-32},{8,32}},
          rotation=90,
          origin={0,68})));

  equation

    connect(inlets,airFlow.inlets);
    connect(outlets,airFlow.outlets);
    connect(metalWall.heatPorts_b,airFlow.heatPorts);

    Q_flow_1=sum(refFlow.heatTransfer.Q_flows);
    Q_flow_2=sum(airFlow.heatPorts.Q_flow);

    Ta_out_ave=sum(airFlow.T_b)/Ncell;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,-60},{100,60}},
          lineColor={28,108,200},
          lineThickness=0.5),
        Line(
          points={{-100,0},{100,0}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{80,-60},{42,-20},{80,20},{40,60}},
          color={255,0,0},
          thickness=1),
        Line(
          points={{24,-60},{-16,-20},{22,20},{-16,60}},
          color={255,0,0},
          thickness=1),
        Line(
          points={{-34,-60},{-72,-20},{-34,20},{-74,60}},
          color={255,0,0},
          thickness=1),                   Text(
            extent={{-58,80},{56,106}},
            lineColor={0,0,127},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="%name")}),                                 Diagram(coordinateSystem(preserveAspectRatio=false)));
  end FinTubeHX;

  package BaseClasses
    extends Modelica.Icons.BasesPackage;

    partial model DryAirCoil "Dry coil model of air side sensible heat transfer for cross flow HX"

      outer TransientVCC.Component.System system;

      import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
      import Modelica.Fluid.Utilities.regStep;

      output SI.HeatFlowRate Q_flows[n] "Total heat transfer rate";

      replaceable package Medium=Modelica.Media.Air.MoistAir;

      parameter Integer n=1 "Number of control volumes";

      parameter Real L=1 "Inertance";

    /*************** Connectors ***************/
      ThermofluidStream.Interfaces.Inlet inlets[n](redeclare each package Medium = Medium,
      each m_flow(min=0))
      "Inlet ports, flow dirction fixed";
      ThermofluidStream.Interfaces.Outlet outlets[n](redeclare each package Medium=Medium,
      each m_flow(max=0))
      "Outlet ports, flow direction fixed";
      TransientVCC.Interfaces.HeatPorts_a heatPorts[n]
      "Heat transfer with metal wall";

    /********************* Fin coil geometry ******************/
      parameter Real Eta_fin_overall=1 "Overall fin efficiency";
      parameter SI.Diameter diameters[n]=ones(n) "Tube outter diameter";
      parameter SI.Length length[n]=ones(n) "Flow path length";
      parameter SI.Area As[n] "Heat transfer area";

      input SI.Area Ac[n](min=1e-10)
        "Effective cross-sectional area (with frost)";

    /********************* Initial conditions ******************/
      parameter Medium.MassFlowRate m_flows_init[n]=ones(n) "Initial air flow rate";
      parameter Medium.Temperature T_out_init[n]=fill(Medium.T_default,n);
      parameter Medium.MassFraction w_out_init[n]=fill(0.01,n);

    /*************** Variables ***************/
      Medium.ThermodynamicState states_a[n] "states of inlet ports";
      Medium.ThermodynamicState states_b[n] "states of outlet ports";
      Medium.Temperature T_a[n](displayUnit="degC")=Medium.temperature(states_a);
      output SI.Temperature T_b[n](displayUnit="degC",start=T_out_init) "Air exit temperatures";
      Medium.MassFraction w_a[n] "Air inlet humidity per unit mass of dry air";
      Medium.MassFraction w_b[n](start=w_out_init) "Air outlet humidity per unit mass of dry air";
      SI.Velocity v[n](each min=0) "Air flow velocity";

    protected
      parameter SI.Temperature T_min=system.T_min;
      parameter SI.Temperature T_max=system.T_max;
      SI.Temperature T_s[n] "Surface temperatures";
      SI.Temperature T_b_forced[n](displayUnit="degC",start=T_out_init) "Air exit temperature under forced convection";
      SI.Temperature T_b_free[n](displayUnit="degC",start=T_out_init) "Under free convection";
      SI.MassFlowRate m_flows[n](start=m_flows_init,each min=0);
      SI.SpecificHeatCapacity cp[n]=Medium.specificHeatCapacityCp(states_a);
      Medium.Density rho[n]=Medium.density(states_a);
      Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states_a);
      Medium.SaturationProperties sat[n];
      Real Ntu[n];
      Real Ntu_free[n];
      SI.HeatFlowRate Q_flows_sen[n] "Sensible heat transfer rates";
      SI.HeatFlowRate Q_flows_free[n] "Heat transfer rate under free convection";
      SI.HeatFlowRate Q_flows_forced[n] "Heat transfer rate under forced convection";
      Medium.AbsolutePressure p_out[n](each start=Medium.p_default);

    /******************** Heat Transfer *************************/

    public
       replaceable model HeatTransfer =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_ConstCoefficient
       constrainedby TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_ConstCoefficient;

       // Forced convection
       HeatTransfer heatTransfer(
        redeclare final package Medium=Medium,
        final n=n,
        final states=states_a,
        final A=As,
        final Ac=Ac,
        final dimension=diameters,
        final length=length,
        final v=v);

       // Free convection
       replaceable model FreeConvection =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_FreeConvection;

       FreeConvection freeConvection(
       redeclare final package Medium=Medium,
        final n=n,
        final states=states_a,
        final A=As,
        final Ac=Ac,
        final dimension=diameters,
        final length=length,
        final v=v,
        final Ts=T_s);

    /******************** Pressure drop *************************/

      replaceable model FrictionalDP =
          TransientVCC.Component.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor;

      FrictionalDP frictionalDP(
      redeclare final package Medium=Medium,
      final m=n,
      final rho=rho,
      final mu=mu,
      final sat=sat,
      final v=v,
      final L=length,
      final dimension=diameters,
      final As=Eta_fin_overall*As,
      final Ac=Ac);

    equation

    /******************* Properties calculation **********************/
      //states_a=Medium.setState_phX(ports_a.p,inStream(ports_a.h_outflow),inStream(ports_a.Xi_outflow));
      states_b={Medium.setState_pTX(p_out[i],T_b[i],{max(0,w_b[i]/(1+w_b[i]))}) for i in 1:n};
      w_a={states_a[i].X[Medium.Water]/states_a[i].X[Medium.Air] for i in 1:n};
      for i in 1:n loop
        sat[i].Tsat=T_a[i];
        sat[i].psat=0.0046;//Medium.saturationPressure(T_a[i]);
      end for;

    /****************** Heat Transfer ******************/
      for i in 1:n loop
        v[i]=m_flows[i]/(rho[i]*Ac[i]);
        Ntu[i]=heatTransfer.alphas[i]*Eta_fin_overall*As[i]/(max(1e-9,m_flows[i])*cp[i]);
        Ntu_free[i]=1*As[i]/(m_flows[i]*cp[i]);
        Q_flows_sen[i]=m_flows[i]*cp[i]*(1-exp(-Ntu[i]))*(T_s[i]-T_a[i]);
        T_b_forced[i]=T_s[i]-(T_s[i]-T_a[i])*exp(-Ntu[i]);
        Q_flows_free[i]=m_flows[i]*cp[i]*(T_b_free[i]-T_a[i]);
        Q_flows_free[i]=1*As[i]*(T_s[i]-(T_b_free[i]+T_a[i])/2);
        Q_flows[i]=regStep(freeConvection.fluidMotion[i]-0.1,Q_flows_free[i],Q_flows_forced[i],1e-4);
        //Q_flows[i]=smooth(0,noEvent(if log10(freeConvection.fluidMotion[i])>1 then Q_flows_free[i] elseif log10(freeConvection.fluidMotion[i])<-2 then Q_flows_forced[i] else (1-log10(freeConvection.fluidMotion[i]))*Q_flows_forced[i]+(log10(freeConvection.fluidMotion[i])+2)*Q_flows_free[i]));
        //T_b[i]=regStep(freeConvection.fluidMotion[i]-0.1,T_b_free[i],T_b_forced[i],1e-4);
        T_b[i]=T_b_forced[i];
      end for;

    /*************** Mass balance *********************/
      inlets.m_flow+outlets.m_flow=zeros(n);

      for i in 1:n loop
        -L*der(m_flows[i])=outlets[i].r-inlets[i].r;
      end for;

    /*************** pressure drop *********************/
      Medium.pressure(states_a) - p_out = frictionalDP.dp;

    /********************* Boundary Conditions ************************/
      T_s=heatPorts.T;
      heatPorts.Q_flow=Q_flows;
      m_flows=inlets.m_flow;
      states_a=inlets.state;
      states_b=outlets.state;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StartTime=225,
          StopTime=2028,
          Tolerance=0.01,
          __Dymola_Algorithm="Dassl"));
    end DryAirCoil;
  end BaseClasses;

  model FinTubeHX_FrostDefrost "Switching between frost formation and melting"

    input Real frostmode;

    input SI.Temperature T_amb;

    // extending refrigerant flow models and tube wall models
    extends DynamicVCC.Components.Units.HX.BaseClasses.PartialHX(final C_metalWall=C_FinTube);

    // Tube and fin geometry
    parameter SI.Mass M_metalWall "Tube wall mass";
    parameter SI.SpecificHeatCapacity cp_metalWall "Tube material specific heat";
    parameter SI.Area Ac_2 "Air side cross-sectional area";
    parameter Real Eta_fin_overall=1 "Overall fin efficiency";
    parameter SI.Area As_2 "Air side heat transfer area";
    parameter SI.Diameter diameter_2;
    parameter SI.Length L_2;
    parameter SI.Mass M_fin "Fin mass";
    parameter SI.SpecificHeatCapacity cp_fin=1 "Fin material specific heat";
    parameter SI.LewisNumber Le=system.Le "Lewis number";

    parameter SI.Height H_hx=1;
    parameter SI.Length L_tube "Single tube length";
    parameter SI.Thickness t_fin=1;
    parameter Real fin_meter=1;
    parameter Integer N_t_prow=1 "# of tubes per row/bank";

    // Initial conditions
    parameter SI.Thickness x_f_init=1e-5;
    parameter SI.Density rho_f_init=30;
    parameter SI.Time T_sample_frost=10 "Sampling period for frost dynamics";

  protected
    parameter SI.HeatCapacity C_FinTube=M_metalWall*cp_metalWall+M_fin*cp_fin "Heat capacity of tube and fin assuming uniform temperature";

    //Initial conditions
  public
    parameter SI.MassFlowRate m_flows_a_init[Ncell]=ones(Ncell) "Initial air flow rate";
    parameter SI.ThermodynamicTemperature Ta_init[Ncell]=fill(Medium_2.reference_T,Ncell) "Initial air exit temperature";
    parameter SI.MassFraction w_a_init[Ncell]=fill(1e-3,Ncell) "Initial air exit humidity ratio";
    SI.ThermodynamicTemperature Ta_out_ave(start=Ta_init[1]);

    /*  Air side  */

    replaceable model HeatTransfer_2 =
        TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_ConstCoefficient
    "Air side heat transfer model";

    replaceable model Friction_2 =
        TransientVCC.Component.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor
    "Air side frictional pressure drop model";

    // Air flow
    replaceable model AirFlow =
        DynamicVCC.Components.Units.SteadyMassFlowPressure.MoistAirCrossFlow;

    AirFlow airFlow(redeclare final package Medium=Medium_2,
    final Ncell=Ncell,
    redeclare final model HeatTransfer=HeatTransfer_2,
    redeclare final model FrictionalDP=Friction_2,
    final Eta_fin_overall=Eta_fin_overall,
    final Le=Le,
    final diameters=fill(diameter_2,Ncell),
    final length=fill(L_2,Ncell),
    final As=fill(As_2/Ncell,Ncell),
    final Ac=Ac_eff,
    final m_flows_init=m_flows_a_init,
    final T_out_init=Ta_init,
    final w_out_init=w_a_init);

    // Air flow connectors
    ThermofluidStream.Interfaces.Inlet inlets[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
          extent={{-8,-32},{8,32}},
          rotation=90,
          origin={0,-68})));

    ThermofluidStream.Interfaces.Outlet outlets[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
          extent={{-8,-32},{8,32}},
          rotation=90,
          origin={0,68})));

    // Frost melting

    replaceable model FrostMelt =
        DynamicVCC.Components.Units.HX.FrostDefrostHX.FrostMelt_Fuzzy_new;

    FrostMelt frostMelt[Ncell](
      each As=As_2/Ncell,
      final airState=airFlow.states_a,
      final T_f_init=Tt_init,
      final tau_f=x_f,
      final rho_f=rho_f,
      final T_t=metalWall.heatPorts_a.T,
      final T_amb=T_amb);

    replaceable model FrostFormation =
        DynamicVCC.Components.Units.HX.FrostDefrostHX.FrostGrowth_Lee;

    FrostFormation frostFormation[Ncell](
    each As=As_2/Ncell,
    final x_f=x_f,
    final rho_f=rho_f,
    final m_flow_dehumid=airFlow.m_flows_dehumid,
    final T_t=metalWall.heatPorts_a.T);

    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatSink_metalWall[Ncell];

  protected
    SI.Area Ac_eff[Ncell] "Varying air cross-sectional area due to frost formation";
    SI.Thickness x_f[Ncell](each start=x_f_init);
    SI.Density rho_f[Ncell](each start=rho_f_init);
    SI.Thickness x_f_formation[Ncell](each start=x_f_init);
    SI.Density rho_f_formation[Ncell](each start=rho_f_init);
  equation

    // Frost dynamics
    when sample(10,T_sample_frost) then
      for i in 1:Ncell loop
        x_f_formation[i]=pre(x_f_formation[i])+frostmode*(frostFormation[i].m_x/pre(rho_f_formation[i])*T_sample_frost)+(if frostmode<1e-8 then 1.0 else 0)*((x_f[i]-x_f_formation[i])/1e-3*T_sample_frost);
        rho_f_formation[i]=pre(rho_f_formation[i])+frostmode*(frostFormation[i].m_rho/pre(x_f_formation[i])*T_sample_frost)+(if frostmode<1e-8 then 1.0 else 0)*((rho_f[i]-rho_f_formation[i])/1e-3*T_sample_frost);

        //x_f_formation[i]=pre(x_f_formation[i])+frostmode*(frostFormation[i].m_x/pre(rho_f_formation[i])*T_sample_frost);
        //rho_f_formation[i]=pre(rho_f_formation[i])+frostmode*(frostFormation[i].m_rho/pre(x_f_formation[i])*T_sample_frost);
        //x_f_formation[i]=(x_f[i]-pre(x_f_formation[i]))/1e-3*T_sample_frost;
        //rho_f_formation[i]=(rho_f[i]-pre(rho_f_formation[i]))/1e-3*T_sample_frost;

      end for;
    end when;

    der(x_f)=frostmode*(x_f_formation-x_f)/1e-4+(1.0-frostmode)*frostMelt.xdot_fuzzy[2];
    der(rho_f)=frostmode*(rho_f_formation-rho_f)/1e-4+(1.0-frostmode)*frostMelt.xdot_fuzzy[7];
    //der(x_f)=frostmode*frostMelt.xdot_fuzzy[2];
    //der(rho_f)=frostmode*frostMelt.xdot_fuzzy[7];

    heatSink_metalWall.Q_flow=Modelica.Fluid.Utilities.regStep(frostmode-0.5,-frostFormation.Q_flow,-frostMelt.Q_flow,1e-3);
    //-frostmode*frostFormation.Q_flow-(1.0-frostmode)*frostMelt.Q_flow;

    connect(inlets,airFlow.inlets);
    connect(outlets,airFlow.outlets);
    connect(metalWall.heatPorts_b,heatSink_metalWall.port);
    //connect(metalWall.heatPorts_b,frostFormation.heatPort_a);
    connect(frostFormation.heatPort,airFlow.heatPorts);
    //connect(metalWall.heatPorts_b,frostMelt.heatPort_a);
    //connect(frostMelt.heatPort_b,airFlow.heatPorts);
    /*
  frostMelt.heatPort_a.T=metalWall.heatPorts_b.T;
  frostFormation.heatPort_a.T=metalWall.heatPorts_b.T;
  metalWall.heatPorts_b.Q_flow+(frostmode*frostFormation.heatPort_a.Q_flow+(1-frostmode)*frostMelt.heatPort_a.Q_flow)=zeros(Ncell);
  airFlow.heatPorts.T=frostFormation.heatPort_b.T;
  airFlow.heatPorts.T=frostMelt.heatPort_b.T;
  airFlow.heatPorts.Q_flow+(frostmode*frostFormation.heatPort_b.Q_flow+(1-frostmode)*frostMelt.heatPort_b.Q_flow)=zeros(Ncell);
*/

    Q_flow_1=sum(refFlow.heatTransfer.Q_flows);
    Q_flow_2=sum(airFlow.heatPorts.Q_flow);

    Ta_out_ave=sum(airFlow.T_b)/Ncell;

    Ac_eff={H_hx*L_tube/Ncell - (t_fin + 2*x_f[i])*L_tube/Ncell*fin_meter*(H_hx - (diameter_2 + 2*x_f[i])*N_t_prow) -
      N_t_prow*(diameter_2 + 2*x_f[i])*L_tube/Ncell for i in 1:Ncell};

  initial equation
     Ac_eff={H_hx*L_tube/Ncell - (t_fin + 2*x_f_init)*L_tube/Ncell*fin_meter*(H_hx - (diameter_2 + 2*x_f_init)*N_t_prow) -
      N_t_prow*(diameter_2 + 2*x_f_init)*L_tube/Ncell for i in 1:Ncell};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,-60},{100,60}},
          lineColor={28,108,200},
          lineThickness=0.5),
        Line(
          points={{-100,0},{100,0}},
          color={0,0,0},
          thickness=1),
        Line(
          points={{80,-60},{42,-20},{80,20},{40,60}},
          color={255,0,0},
          thickness=1),
        Line(
          points={{24,-60},{-16,-20},{22,20},{-16,60}},
          color={255,0,0},
          thickness=1),
        Line(
          points={{-34,-60},{-72,-20},{-34,20},{-74,60}},
          color={255,0,0},
          thickness=1),                   Text(
            extent={{-58,80},{56,106}},
            lineColor={0,0,127},
            lineThickness=0.5,
            fillColor={0,0,255},
            fillPattern=FillPattern.Solid,
            textString="%name")}),                                 Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StartTime=225,
        StopTime=2320,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end FinTubeHX_FrostDefrost;

  package Tests
    extends Modelica.Icons.ExamplesPackage;
    model Test_Cycle_FrostingDefrostingOD28 "Cycling of frosting and defrsoting for an entire experiment test, validation dataset 01282022"
      extends Modelica.Icons.Example;

       inner TransientVCC.Component.System system(
      p_max=45e5,
      p_min=2e5,
      h_max=4.7e5,
      h_min=1.11e5,
      T_max=340,
      T_min=250,
      EnableReverseFlow=EnableReverseFlow);

      parameter Integer Ncell=30;

      //replaceable package Medium_CP=TransientVCC.Media.CoolProp.R410a;
      replaceable package Medium_NN=TransientVCC.Media.R410A_ph;
      package Medium_1=Medium_NN;
      package Medium_2=Modelica.Media.Air.MoistAir;

    /***************************** HX Geometry *******************************/

      DynamicVCC.Examples.Tests.FinTubeGeo Geo_OD(
        D_o=0.0074,
        D_i=0.0068,
        cp_tube=385,
        rho_tube=8900,
        cp_fin=900,
        rho_fin=2700,
        L_tube=2.7026,
        N_tuberow=2,
        N_t_prow=44,
        N_circuits=8,
        pf=0.0012,
        pt=0.0216,
        pl=0.0187,
        fin_meter=787.4016,
        t_fin=9.9060e-05);

      DynamicVCC.Examples.Tests.FinTubeGeo Geo_ID(
        D_o=0.01,
        D_i=0.0086,
        cp_tube=900,
        rho_tube=2700,
        cp_fin=900,
        rho_fin=2700,
        L_tube=0.4521,
        N_tuberow=3,
        N_t_prow=28,
        N_circuits=6,
        pf=0.0016,
        pt=0.0254,
        pl=0.0191,
        fin_meter=570.8661,
        t_fin=1.0668e-04);

    /***************************** Constants *******************************/

      parameter SI.AbsolutePressure p_scale=1e5;
      parameter SI.SpecificEnthalpy h_scale=1e5;
      parameter SI.Temperature T_scale=273;
      parameter SI.CoefficientOfHeatTransfer alpha_nominal_ID=10000;
      parameter SI.CoefficientOfHeatTransfer alpha_nominal_OD=5000;
      parameter SI.LewisNumber Le=0.854 "Lewis number";
      parameter Integer N_pipe=5; //Discharge pipe line
      parameter SI.Time samplePeriod_exv=90;
      parameter SI.Time startTime=0 "Simulation start time";

    /***************************** Initialization *******************************/

      parameter Real p_init_ID[Ncell]=fill(7.68e5,Ncell);
      parameter Real h_init_ID[Ncell]=fill(4.41e5,Ncell);
      parameter Real Tt_init_ID[Ncell]=fill(291.4,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_ID[Ncell]=fill(1e-3/Ncell,Ncell);
      parameter SI.Temperature T_a_init_ID[Ncell]=fill(292,Ncell);
      parameter SI.MassFraction w_a_init_ID[Ncell]=fill(0.004,Ncell);

      parameter Real p_init_OD[Ncell]=fill(7.6e5,Ncell);
      parameter Real h_init_OD[Ncell]=fill(4.22e5,Ncell);
      parameter Real Tt_init_OD[Ncell]=fill(272,Ncell);
      parameter SI.MassFlowRate m_dot_air_ini_OD[Ncell]=fill(1e-4/Ncell,Ncell);
      parameter SI.Temperature T_a_init_OD[Ncell]=fill(272,Ncell);
      parameter SI.MassFraction w_a_init_OD[Ncell]=fill(0.0023,Ncell);

      import TransientVCC.Component.Types.ModelStructure;
      parameter Boolean SteadyState_init=false;
      parameter Boolean EnableReverseFlow=true;
      parameter Boolean useLumpedPressure=false;
      parameter Boolean use_I_flows=if (useLumpedPressure and modelStructure==ModelStructure.av_vb) then false else true;
      parameter SI.MassFlowRate m_flows_init_1[Ncell+1]=fill(1e-3,Ncell+1);
      parameter TransientVCC.Component.Types.ModelStructure modelStructure=TransientVCC.Component.Types.ModelStructure.av_vb;

      parameter SI.Frequency speed_init=1;
      parameter SI.Power Pwr_init=1;
      parameter Real opening_init=0.1;
      parameter SI.Pressure dp_init=23e5;

      parameter SI.Thickness x_f_init=1e-5;
      parameter SI.Density rho_f_init=30;

    /***************************** Heat Transfer and Pressure Drop *******************************/

    /*
  replaceable model Condensation =TransientVCC.Component.Pipes_newcell.BaseClasses.HeatTransfer.Correlations.Condensation_Shah (
  pc=4901200,
  alpha0=alpha_nominal_ID);
*/

      replaceable model Condensation =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
      alpha0=alpha_nominal_ID);

      replaceable model Evaporation =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
      alpha0=alpha_nominal_OD);

      replaceable model Singlephase =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.Correlations.SinglePhase_Gnielinski (
      alpha0=alpha_nominal_OD);

     /*
  replaceable model Singlephase =TransientVCC.Component.Pipes_newcell.BaseClasses.HeatTransfer.Correlations.Constant (
  alpha0=3000);
  */
      replaceable model HeatTransfer_1_OD =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.HeatTransferPhaseZones (
      redeclare model LiquidZone=Singlephase,
      redeclare model VaporZone=Singlephase,
      redeclare model TwoPhaseZone=Evaporation);

     replaceable model HeatTransfer_2_OD =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_Wang (
      Nrow=Geo_OD.N_tuberow,
      pf=Geo_OD.pf,
      pl=Geo_OD.pl,
      pt=Geo_OD.pt,
      t_fin=Geo_OD.t_fin,
      final alpha0=45);

    /*
  replaceable model HeatTransfer_2_OD=TransientVCC.Component.Pipes_newcell.BaseClasses.HeatTransfer.AirCoilHeatTransfer_ConstCoefficient (
  final alpha0=50);
*/
     replaceable model Friction_1_OD =
          TransientVCC.Component.Pipes.BaseClasses.Friction.Correlations.Constant (                                  f0=0.1);

    /*
 replaceable model Friction_2_OD=TransientVCC.Component.Pipes_newcell.BaseClasses.Friction.AirCoilDP_Wang (
      Nt=Geo_OD.N_tuberow,
      D=Geo_OD.D_o,
      pf=Geo_OD.pf,
      pl=Geo_OD.pl,
      pt=Geo_OD.pt,
      t_fin=Geo_OD.t_fin,
      final dp_nominal=10);
*/

     replaceable model Friction_2_OD =
          TransientVCC.Component.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
     f0=0.12)
     "Outdoor air side pressure drop";

    /*
replaceable model HeatTransfer_1_ID =TransientVCC.Component.Pipes_newcell.BaseClasses.HeatTransfer.ConstHTC (
alpha0=fill(alpha_nominal_ID,Ncell));
*/

      replaceable model HeatTransfer_1_ID =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.HeatTransferPhaseZones (
      redeclare model LiquidZone=Singlephase,
      redeclare model VaporZone=Singlephase,
      redeclare model TwoPhaseZone=Condensation);

     //replaceable model Friction_1_ID=TransientVCC.Component.Pipes_newcell.BaseClasses.Friction.Correlations.TwoPhase_SinglePhase;
     replaceable model Friction_1_ID =
          TransientVCC.Component.Pipes.BaseClasses.Friction.Correlations.Constant (                                  f0=0.1);

    /*

 replaceable model HeatTransfer_1_ID =
  TransientVCC.Component.Units.HeatTransfer.HeatTransferCorrelation.ConstTwoPhaseGnielinskiDittusBoelter (
  h0=HTC_nominal_ID,h_2phase=HTC_nominal_ID)
 "Indoor refrigerant side heat transfer model";
 */

     replaceable model HeatTransfer_2_ID =
          TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.AirCoilHeatTransfer_Wang (
      Nrow=Geo_ID.N_tuberow,
      pf=Geo_ID.pf,
      pl=Geo_ID.pl,
      pt=Geo_ID.pt,
      t_fin=Geo_ID.t_fin,
      final alpha0=45)
      "Indoor air side heat transfer model";

    /*
  replaceable model Friction_2_ID=TransientVCC.Component.Pipes_newcell.BaseClasses.Friction.AirCoilDP_Wang (
      Nt=Geo_ID.N_tuberow,
      D=Geo_ID.D_o,
      pf=Geo_ID.pf,
      pl=Geo_ID.pl,
      pt=Geo_ID.pt,
      t_fin=Geo_ID.t_fin,
      final dp_nominal=20);
*/

     replaceable model Friction_2_ID =
          TransientVCC.Component.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor (
     f0=0.1)
     "Indoor air side pressure drop";

    /********************************* Control ******************************************/

      Modelica.Blocks.Continuous.FirstOrder Comp_speed(T=5,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=speed_init);
      Modelica.Blocks.Continuous.FirstOrder EXV_opening(T=3,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=opening_init);
      Modelica.Blocks.Continuous.FirstOrder RV_opening(T=3,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=1.0);
      Modelica.Blocks.Continuous.FirstOrder fan_OD_speed(T=10,initType=Modelica.Blocks.Types.Init.InitialOutput,y_start=1);

      Modelica.Blocks.Nonlinear.Limiter Comp_limiter(uMax=110,uMin=0.85);
      Modelica.Blocks.Nonlinear.Limiter exv_limiter(uMax=1.0,uMin=0.01);

    /*
  Modelica.Blocks.Discrete.Sampler sampler_sh(samplePeriod=samplePeriod_exv,startTime=startTime);
  Modelica.Blocks.Discrete.ZeroOrderHold ZOH_exv(samplePeriod=samplePeriod_exv,startTime=startTime);

  Modelica.Blocks.Continuous.PI exv_controller(
  k=-4e-3,
  T=100,
  y_start=0.1,
  x_start=-25,
  initType=Modelica.Blocks.Types.Init.InitialOutput);

  Modelica.Blocks.Nonlinear.Limiter limiter_exv(uMax=1,uMin=0.1);
  */

      //Modelica.Blocks.Discrete.Sampler sampler_comp(samplePeriod=10,startTime=startTime);
      //Modelica.Blocks.Discrete.ZeroOrderHold ZOH_comp(samplePeriod=10,startTime=startTime);
    /*
  TransientVCC.Simulations.HeatPump.CompressorPIController compressor_controller(
  samplePeriod=5,
  startTime=startTime);
     */

        Modelica.Blocks.Sources.CombiTimeTable BC_ID(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_ID",fileName="C:/Jiacheng Ma/CarrierGreenspeed/BC_ID.mat",columns=2:5);
        Modelica.Blocks.Sources.CombiTimeTable BC_OD(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_OD",fileName="C:/Jiacheng Ma/CarrierGreenspeed/BC_OD.mat",columns=2:5);
        Modelica.Blocks.Sources.CombiTimeTable Mea_ID(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea_ID",fileName="C:/Jiacheng Ma/CarrierGreenspeed/Mea_ID.mat",columns=2:6);
        Modelica.Blocks.Sources.CombiTimeTable Mea_OD(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea_OD",fileName="C:/Jiacheng Ma/CarrierGreenspeed/Mea_OD.mat",columns=2:6);
        Modelica.Blocks.Sources.CombiTimeTable Comp(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Comp",fileName="C:/Jiacheng Ma/CarrierGreenspeed/Comp.mat",columns=2:15);
        Modelica.Blocks.Sources.CombiTimeTable EXVopen(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="data",fileName="C:/Jiacheng Ma/CarrierGreenspeed/EXVopen.mat",columns=2:2);

    /***************************** Components *******************************/

      DynamicVCC.Components.Units.MassFlowDevices.Compressor.Map10Coefficient compressor(
        redeclare final package Medium = Medium_1,
        Vs=2.3765e-05,
        T_amb=T_amb,
        m_flow_init=m_flows_init_1[1],
        h_dis_init=h_init_ID[1],
        h_suc_init=h_init_OD[Ncell],
        p_dis_init=p_init_ID[1],
        p_suc_init=p_init_OD[Ncell],
        speed_init=speed_init,
        Pwr_init=Pwr_init) annotation (Placement(transformation(
            extent={{-26,-26},{26,26}},
            rotation=90,
            origin={148,18})));

    /*
  TransientVCC.Component.Units.MassFlowDevices_Explicit.Compressor compressor(
  redeclare final package Medium=Medium_1,
  final Vs=2.3765e-5,
  T_amb=T_amb,
  L=1,
  m_flow_init=m_flows_init_1[1],
  speed_init=speed_init,
  SteadyState_init=true);
  */

      MassFlowDevices.Valve.ElectronicExpansionValve exv(
        redeclare package Medium = Medium_1,
        final Av=2.5447e-06,
        EnableReverseFlow=EnableReverseFlow,
        p_a_init=p_init_ID[Ncell],
        p_b_init=p_init_OD[1],
        opening_init=opening_init,
        dp_nominal=23e5,
        dp_init=dp_init,
        m_flow_init=m_flows_init_1[1]) annotation (Placement(transformation(
            extent={{14,-14},{-14,14}},
            rotation=90,
            origin={-134,22})));

      HX.FinTubeHX IndoorCoil(
        redeclare final package Medium_1 = Medium_1,
        redeclare final package Medium_2 = Medium_2,
        redeclare final model HeatTransfer_1 = HeatTransfer_1_ID,
        redeclare final model HeatTransfer_2 = HeatTransfer_2_ID,
        redeclare final model Friction_1 = Friction_1_ID,
        redeclare final model Friction_2 = Friction_2_ID,
        final modelStructure=modelStructure,
        final EnableReverseFlow=EnableReverseFlow,
        final useLumpedPressure=useLumpedPressure,
        final use_I_flows=use_I_flows,
        final Ncell=Ncell,
        As_1=Geo_ID.HTA_r,
        Ac_1=Geo_ID.Ac_r,
        L_1=Geo_ID.L_circuit,
        M_metalWall=Geo_ID.M_tube,
        cp_metalWall=Geo_ID.cp_tube,
        diameter_2=Geo_ID.D_o,
        diameter_1=Geo_ID.D_i,
        L_2=Geo_ID.pl,
        Eta_fin_overall=Geo_ID.Eta_fin_overall,
        Le=Le,
        As_2=Geo_ID.HTA_e,
        Ac_2=Geo_ID.Ac_e,
        M_fin=Geo_ID.M_fin,
        cp_fin=Geo_ID.cp_fin,
        p_init=p_init_ID,
        h_init=h_init_ID,
        Tt_init=Tt_init_ID,
        m_flows_init=m_flows_init_1,
        p_scale=p_scale,
        h_scale=h_scale,
        T_scale=T_scale,
        w_a_init=w_a_init_ID,
        SteadyState_init=SteadyState_init,
        m_flows_a_init=m_dot_air_ini_ID,
        Ta_init=T_a_init_ID) annotation (Placement(transformation(extent={{16,30},{-56,102}})));
    /*
    replaceable model FrostFormation=TransientVCC.Component.HeatFlow.FrostGrowthModel.FrostGrowth_Lee (
    final x_f_init=x_f_init,
    final rho_f_init=rho_f_init);
*/
      Real frostmode;
      DynamicVCC.Components.Units.SteadyMassFlowPressure.FinTubeHX_FrostDefrost OutdoorCoil(
        redeclare final package Medium_1 = Medium_1,
        redeclare final package Medium_2 = Medium_2,
        redeclare final model HeatTransfer_1 = HeatTransfer_1_OD,
        redeclare final model HeatTransfer_2 = HeatTransfer_2_OD,
        redeclare final model Friction_1 = Friction_1_OD,
        redeclare final model Friction_2 = Friction_2_OD,
        final modelStructure=modelStructure,
        final EnableReverseFlow=EnableReverseFlow,
        final useLumpedPressure=useLumpedPressure,
        final use_I_flows=use_I_flows,
        final Ncell=Ncell,
        As_1=Geo_OD.HTA_r,
        Ac_1=Geo_OD.Ac_r,
        L_1=Geo_OD.L_circuit,
        M_metalWall=Geo_OD.M_tube,
        cp_metalWall=Geo_OD.cp_tube,
        diameter_2=Geo_OD.D_o,
        diameter_1=Geo_OD.D_i,
        L_2=Geo_OD.pl,
        Eta_fin_overall=Geo_OD.Eta_fin_overall,
        Le=Le,
        As_2=Geo_OD.HTA_e,
        Ac_2=Geo_OD.Ac_e,
        M_fin=Geo_OD.M_fin,
        cp_fin=Geo_OD.cp_fin,
        p_init=p_init_OD,
        h_init=h_init_OD,
        Tt_init=Tt_init_OD,
        m_flows_init=m_flows_init_1,
        p_scale=p_scale,
        h_scale=h_scale,
        T_scale=T_scale,
        w_a_init=w_a_init_OD,
        SteadyState_init=SteadyState_init,
        m_flows_a_init=m_dot_air_ini_OD,
        Ta_init=T_a_init_OD,
        final L_tube=Geo_OD.L_tube,
        final H_hx=Geo_OD.H_hx,
        final t_fin=Geo_OD.t_fin,
        final fin_meter=Geo_OD.fin_meter,
        final N_t_prow=Geo_OD.N_t_prow,
        final frostmode=frostmode,
        final T_amb=T_amb) annotation (Placement(transformation(extent={{-52,-92},{20,-20}})));

      MassFlowDevices.Valve.ReversingValve reversingValve(
        redeclare package Medium = Medium_1,
        m_flow_init=m_flows_init_1[1],
        p_dis_init=p_init_ID[1],
        p_suc_init=p_init_OD[Ncell],
        A_open=1.0e-4,
        dp_nominal=0.2e5,
        C_Hd=1,
        C_Hs=0.7,
        tau=1) annotation (Placement(transformation(
            extent={{-19,-19},{19,19}},
            rotation=-90,
            origin={79,-7})));

    /*
  TransientVCC.Component.Units.MassFlowDevices_Explicit.ReversingValve reversingValve(
    redeclare package Medium=Medium_1,
    final SteadyState_init=false,
    tau=1,
    opening_min=1e-5,
    A_open=8.77e-4,
    m_flow_init=m_flows_init_1[1],
    L=1);
*/
      MassFlowDevices.Accumulator accumulator(
        redeclare package Medium = Medium_1,
        m_flow_init=m_flows_init_1[1],
        V=0.0058,
        V_f_init=0.5*0.0058,
        p_init=p_init_OD[Ncell],
        h_init=h_init_OD[Ncell],
        SteadyState_init=SteadyState_init) annotation (Placement(transformation(extent={{108,-78},{142,-44}})));

      SI.MassFlowRate m_OD_air(start=sum(m_dot_air_ini_OD));
      SI.MassFlowRate m_OD_air_SP;
      SI.Temperature T_amb;

      Modelica.Fluid.Sources.MassFlowSource_T fan_ID[Ncell](
        redeclare each final package Medium = Medium_2,
        each use_m_flow_in=true,
        each use_T_in=true,
        each use_X_in=true,
        each nPorts=1) annotation (Placement(transformation(extent={{-88,22},{-68,42}})));

    /*
  Modelica.Fluid.Sources.MassFlowSource_T fan_OD[Ncell](
    redeclare each final package Medium = Medium_2,
    each use_m_flow_in=true,
    each use_T_in=true,
    each use_X_in=true,
    each nPorts=1) annotation (Placement(transformation(extent={{-98,-92},{-78,-72}})));
*/

      package Medium_stream=ThermofluidStream.Media.myMedia.Air.MoistAir;

      inner ThermofluidStream.DropOfCommons dropOfcommons;

      ThermofluidStream.Boundaries.Source airSource(
      redeclare final package Medium=Medium_stream,
      final p0_par=Medium_2.p_default,
      final temperatureFromInput=true,
      final xiFromInput=true);

      ThermofluidStream.Topology.SplitterN splitter(
      redeclare final package Medium=Medium_stream,
      final N=Ncell,
      final L=20);

      ThermofluidStream.Topology.JunctionN junction(
      redeclare final package Medium=Medium_stream,
      final N=Ncell,
      final L=20,
      final m_flow_eps=1e-8);

      DynamicVCC.Components.Units.SteadyMassFlowPressure.Fan_new fan_OD(
        final Ncell=Ncell,
        V_dot_init=1e-2,
        L=20);

      ThermofluidStream.Boundaries.Sink airSink(
      redeclare final package Medium=Medium_stream,
      final p0_par=Medium_2.p_default,
      final L=20);

      Modelica.Fluid.Sources.Boundary_pT sink_ID(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{-106,88},{-86,108}})));
      //Modelica.Fluid.Sources.Boundary_pT sink_OD(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{-92,-24},{-72,-4}})));

      Sensors.T_Superheat superheat(redeclare package Medium = Medium_1) annotation (Placement(transformation(extent={{40,-86},{60,-66}})));
      Sensors.T_Superheat subcooling(redeclare package Medium = Medium_1) annotation (Placement(transformation(extent={{-148,78},{-128,98}})));

      Sensors.Temperature_grid T_supply_ID(redeclare package Medium = Medium_2, nPorts=Ncell) annotation (Placement(transformation(extent={{20,84},{40,104}})));
      /*
  Component.Units.Sensors.Temperature_grid T_supply_OD(
  redeclare package Medium=Medium_2,
  nPorts=Ncell) annotation (Placement(transformation(extent={{-120,-42},{-100,-22}})));
*/
      TransientVCC.Component.Pipes.ConnectingPipe pipe(
        redeclare final package Medium_1 = Medium_1,
        redeclare final model HeatTransfer_1 =
            TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.ConstHTC (                                   alpha0=100),
        redeclare final model HeatTransfer_2 =
            TransientVCC.Component.Pipes.BaseClasses.HeatTransfer.HeatLossAmbient (                                   alpha0=10),
        redeclare final model Friction_1 = Friction_1_OD,
        final modelStructure=TransientVCC.Component.Types.ModelStructure.av_vb,
        final EnableReverseFlow=false,
        final SteadyState_init=SteadyState_init,
        final useLumpedPressure=useLumpedPressure,
        final use_I_flows=use_I_flows,
        final Ncell=N_pipe,
        final T_amb=T_amb,
        final p_scale=p_scale,
        final h_scale=h_scale,
        final T_scale=T_scale,
        d_i=0.0110744,
        d_o=0.0127,
        L_tube=1.0922,
        p_init=fill(p_init_ID[1], N_pipe),
        h_init=fill(h_init_ID[1], N_pipe),
        Tt_init=fill(Tt_init_ID[1], N_pipe),
        m_flows_init=fill(m_flows_init_1[1], N_pipe + 1)) annotation (Placement(transformation(extent={{130,60},{110,80}})));

    /*

  Modelica.Fluid.Sources.MassFlowSource_h source(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_m_flow_in=true,
  use_h_in=true);

  Modelica.Fluid.Sources.MassFlowSource_h sink(
  redeclare package Medium=Medium_1,
  nPorts=1,
  use_m_flow_in=true);
*/

    equation

    /*
   for i in 1:Ncell loop
     fan_OD[i].m_flow_in=m_OD_air/Ncell;
     fan_OD[i].X_in={0.0022,1-0.0022};
     fan_OD[i].T_in=293;
   end for;
*/

       der(m_OD_air)=(m_OD_air_SP-m_OD_air)/1e-3;
       m_OD_air_SP=max(1e-5,BC_OD.y[2]);

    /*
  for i in 1:Ncell loop
     fan_OD[i].m_flow_in=m_OD_air/Ncell;
     fan_OD[i].X_in={BC_OD.y[3],1-BC_OD.y[3]};
     fan_OD[i].T_in=BC_OD.y[1];
  end for;
*/
      airSource.T0_var=BC_OD.y[1];
      airSource.xi_var={BC_OD.y[3]};

      //fan_OD_speed.u=if time<2457 then Comp.y[13] elseif time<2805 then 1e-3 elseif time<10180 then Comp.y[13] elseif time<10825 then 1e-3 else Comp.y[13];
      fan_OD_speed.u=if time<2457 then max(1,Comp.y[13]) elseif time<2805 then 1 elseif time<10185 then max(1,Comp.y[13]) elseif time<10825 then 1 else max(1,Comp.y[13]);

      connect(fan_OD_speed.y,fan_OD.speed);

      connect(airSource.outlet,fan_OD.inlet);
      connect(fan_OD.outlet,splitter.inlet);
      connect(splitter.outlets,OutdoorCoil.inlets);
      connect(OutdoorCoil.outlets,junction.inlets);
      connect(junction.outlet,airSink.inlet);

      T_amb=BC_OD.y[1];

    /*
  for i in 1:Ncell loop
   fan_ID[i].m_flow_in=1e-3/Ncell;
   fan_ID[i].T_in=294;
   fan_ID[i].X_in={0.0073,1-0.0073};
  end for;
*/

      for i in 1:Ncell loop
         fan_ID[i].m_flow_in=BC_ID.y[2]/Ncell;
         fan_ID[i].T_in=BC_ID.y[1];
         fan_ID[i].X_in={BC_ID.y[3],1-BC_ID.y[3]};
      end for;

       Comp_speed.u=Comp.y[4];
       connect(Comp_speed.y,Comp_limiter.u);
       connect(compressor.speed,Comp_limiter.y);

       RV_opening.u=if time<2519 then 1.0 elseif time <2802 then 0.0 elseif time<10250 then 1.0 elseif time<10820 then 0.0 else 1.0;
       connect(RV_opening.y,reversingValve.opening_SP);

       //EXV_opening.u=if time<617 then 0.1 elseif time<685 then 0.4 elseif time<804 then 0.2 elseif time<2457 then Comp.y[9]-0.01 elseif time<2510 then 0.1+(time-2457)/53*0.9 elseif time<2801 then 1.0 elseif time<2845 then 0.6 elseif time<3000 then 0.26 elseif time<6450 then Comp.y[9]-0.01 elseif time<10000 then Comp.y[9]+(time-6450)/3550*0.01 elseif time<10180 then Comp.y[9] elseif time<10250 then 0.1+(time-10180)/70*0.9 elseif time<10820 then 1.0 elseif time<10840 then 0.6 elseif time<11000 then 0.26 else Comp.y[9];
       //EXV_opening.u=if time<617 then 0.1 elseif time<685 then 0.4 elseif time<804 then 0.2 elseif time<1769 then Comp.y[9]-0.015 elseif time<2350 then 0.093 elseif time<2457 then Comp.y[9]-0.015 elseif time<2510 then 0.1+(time-2457)/53*0.9 elseif time<2801 then 1.0 elseif time<2845 then 0.6 elseif time<3000 then 0.26 elseif time<5000 then Comp.y[9]-0.015 elseif time<5117 then 0.08586 elseif time<5300 then Comp.y[9]-0.015 elseif time<5310 then Comp.y[9]-0.01+(time-5300)*0.0005 elseif time<5580 then Comp.y[9]-0.005 elseif time<5810 then 0.0865 elseif time<6450 then Comp.y[9]-0.005 elseif time<10000 then Comp.y[9]-0.005+(time-6450)/3550*0.01 elseif time<10180 then Comp.y[9] elseif time<10250 then 0.1+(time-10180)/70*0.9 elseif time<10820 then 1.0 elseif time<10840 then 0.6 elseif time<11000 then 0.26 else Comp.y[9];
       //EXV_opening.u=if time<617 then 0.1 elseif time<685 then 0.4 elseif time<804 then 0.2 elseif time<1769 then Comp.y[9]-0.015 elseif time<2350 then 0.093 elseif time<2457 then Comp.y[9]-0.015 elseif time<2510 then 0.1+(time-2457)/53*0.9 elseif time<2801 then 1.0 elseif time<2845 then 0.6 elseif time<3000 then 0.26 elseif time<4618 then Comp.y[9]-0.012 elseif time<4790 then 0.092 elseif time<5000 then Comp.y[9]-0.012 elseif time<5117 then 0.08886 elseif time<5300 then Comp.y[9]-0.012 elseif time<5310 then Comp.y[9]-0.012+(time-5300)*0.0005 elseif time<5580 then Comp.y[9]-0.007 elseif time<6450 then Comp.y[9]-0.005 elseif time<10000 then Comp.y[9]-0.005+(time-6450)/3550*0.01 elseif time<10180 then Comp.y[9] elseif time<10250 then 0.1+(time-10180)/70*0.9 elseif time<10820 then 1.0 elseif time<10840 then 0.6 elseif time<11000 then 0.26 else Comp.y[9];
       EXV_opening.u=EXVopen.y[1];
       connect(EXV_opening.y,exv_limiter.u);
       connect(exv_limiter.y,exv.opening);

       frostmode=smooth(0,noEvent(if time<2440 then 1.0 elseif time<2441 then 1.0-(time-2440) elseif time<2810 then 0.0 elseif time<2811 then (time-2810) elseif time<10160 then 1.0 elseif time<10161 then 1.0-(time-10160) elseif time<10820 then 0.0 elseif time<10821 then time-10820 else 1.0));
    /*
  connect(source.ports[1],pipe.port_a1);
  connect(pipe.port_b1,reversingValve.port_a);
  connect(accumulator.port_b,sink.ports[1]);
  source.m_flow_in=0.0474;
  source.h_in=4.354e5;
  sink.m_flow_in=-0.0474;
  */

      connect(IndoorCoil.port_b1, exv.port_a) annotation (Line(points={{-56,66},{-134,66},{-134,36}}, color={0,127,255}));
      connect(exv.port_b, OutdoorCoil.port_a1) annotation (Line(points={{-134,8},{-134,-56},{-52,-56}},   color={0,127,255}));
      connect(fan_ID.ports[1], IndoorCoil.ports_a2) annotation (Line(points={{-68,32},{-68,12},{-20,12},{-20,41.52}}, color={0,127,255}));
      connect(sink_ID.ports, IndoorCoil.ports_b2) annotation (Line(points={{-86,98},{-62,98},{-62,108},{-20,108},{-20,90.48}}, color={0,127,255}));
      connect(superheat.port, OutdoorCoil.port_b1) annotation (Line(points={{50,-86},{50,-92},{28,-92},{28,-56},{20,-56}}, color={0,127,255}));
      connect(subcooling.port, IndoorCoil.port_b1) annotation (Line(points={{-138,78},{-138,66},{-56,66}}, color={0,127,255}));
      connect(T_supply_ID.ports, IndoorCoil.ports_b2) annotation (Line(points={{30,84.4},{30,76},{50,76},{50,110},{-20,110},{-20,90.48}}, color={0,127,255}));
      //connect(T_supply_OD.ports, OutdoorCoil.ports_b2) annotation (Line(points={{-110,-41.6},{-110,-44},{-62,-44},{-62,-14},{-16,-14},{-16,-31.52}}, color={0,127,255}));
      connect(compressor.port_a, accumulator.port_b) annotation (Line(points={{148,-8},{148,-61},{142,-61}}, color={0,127,255}));
      connect(compressor.port_b, pipe.port_a1) annotation (Line(points={{148,44},{148,70},{130,70}}, color={0,127,255}));
      connect(pipe.port_b1, reversingValve.port_a) annotation (Line(points={{110,70},{102,70},{102,12},{79,12}},                    color={0,127,255}));
      connect(reversingValve.port_b, accumulator.port_a) annotation (Line(points={{79,-26},{78,-26},{78,-61},{108,-61}}, color={0,127,255}));
      connect(reversingValve.port_b2, IndoorCoil.port_a1) annotation (Line(points={{66.46,-26},{66.46,-32},{26,-32},{26,66},{16,66}}, color={0,127,255}));
      connect(OutdoorCoil.port_b1, reversingValve.port_a2) annotation (Line(points={{20,-56},{92,-56},{92,-26},{91.54,-26}}, color={0,127,255}));
      annotation (experiment(
          StartTime=615,
          StopTime=11830,
          Tolerance=0.005,
          __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(extent={{-240,-140},{240,160}})),
        Icon(coordinateSystem(extent={{-240,-140},{240,160}})));
    end Test_Cycle_FrostingDefrostingOD28;
  end Tests;
end SteadyMassFlowPressure;
