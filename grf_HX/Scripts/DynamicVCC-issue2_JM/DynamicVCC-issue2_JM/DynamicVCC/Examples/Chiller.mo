within DynamicVCC.Examples;
package Chiller "90-ton centrifugal chiller"
  extends Modelica.Icons.VariantsPackage;

  model RefFlow1D_UniformPressure "1D flow neglecting pressure drop and eliminate momentum balances (no flow cells)"

    import DynamicVCC.Components.Types.ModelStructure
    "determine number of flow cells (momentum balances) and boundary conditions at port_b";

    //extending connectors at inlet and outlet
    extends DynamicVCC.Interfaces.PartialTwoPort(port_a(m_flow(start=m_flows_init[1])),port_b(m_flow(start=-m_flows_init[Ncell])));

    //Volume cells
    extends DynamicVCC.Components.Pipes.BaseClasses.PartialVolumeCell(final n=Ncell, final V={Ac[i]*L[i] for i in 1:n});

    // Geometry
    parameter SI.Length L[n] "Length of refrigerant flow path";
    parameter SI.Area Ac[n] "Cross flow area";
    parameter SI.Area As[n] "Surface area (heat transfer area)";
    parameter SI.Diameter diameters[n] "Pipe diameter";

    // Discretization
    parameter Integer Ncell(min=1)=2 "Number of volume cells";

    parameter Boolean useLumpedPressure=false
    "=true neglect pressure drop across flow path, use a lumped pressure for all volumes";
    parameter ModelStructure modelStructure=ModelStructure.av_b;
    final parameter Integer nFlowCell=if useLumpedPressure then nFLumped else nFDistributed;
    final parameter Integer nFDistributed=if modelStructure==ModelStructure.a_v_b       then n+2 elseif modelStructure==ModelStructure.av_vb       then n else n+1
    "Number of Flow cells for distributed volumes";
    final parameter Integer nFLumped=if modelStructure==ModelStructure.a_v_b       then 3 else 2
    "Number of Flow cells under lumped pressure";
    final parameter Integer iLumped=integer(n/2)+1;

    //Initialization
    parameter Medium.MassFlowRate m_flows_init[n+1]=ones(n+1);

    // Flow quantities
    Medium.MassFlowRate m_flows[n+1](start=m_flows_init, each min=if EnableReverseFlow then -Modelica.Constants.inf else 0)
      "Mass flow rate of each volume cell boundary";

    Medium.EnthalpyFlowRate H_flows[n+1]
    "Enthalpy flow rate of each volume cell boundary";

    SI.Velocity v_m[n]={ 0.5*(m_flows[i]+m_flows[i+1])/rho[i]/Ac[i] for i in 1:n}
    "Mean velocity of each volume cell";

    // Tube wall heat transfer
    replaceable model HeatTransfer=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.PartialFlowHeatTransfer;

    HeatTransfer heatTransfer(
    redeclare final package Medium=Medium,
    final n=n,
    final states=states,
    final A=As,
    final dimension=diameters,
    final length=L,
    final v=v_m);

    replaceable model FlowCells=DynamicVCC.Components.Pipes.BaseClasses.StaggeredGridMomentum; // Not used

    DynamicVCC.Interfaces.HeatPorts_a heatPorts[n];

  protected
    Medium.ThermodynamicState state_a "State at port_a";
    Medium.ThermodynamicState state_b "State at port_b";
  equation

    // Source terms for mass and energy balances
    Q_flows=heatTransfer.Q_flows;
    m_flows_cell=m_flows[1:n]-m_flows[2:n+1];
    H_flows_cell=H_flows[1:n]-H_flows[2:n+1];

    for i in 2:n loop
      H_flows[i]=semiLinear(m_flows[i], Medium.specificEnthalpy(states[i-1]), Medium.specificEnthalpy(states[i]));
    end for;
    H_flows[1]=homotopy(semiLinear(port_a.m_flow, inStream(port_a.h_outflow), Medium.specificEnthalpy(states[1])),
    port_a.m_flow*inStream(port_a.h_outflow));
    H_flows[n+1]=homotopy(-semiLinear(port_b.m_flow, inStream(port_b.h_outflow), Medium.specificEnthalpy(states[n])),
    -port_b.m_flow*Medium.specificEnthalpy(states[n]));

    // Boundary conditions
    m_flows[1]=port_a.m_flow;
    m_flows[n+1]=-port_b.m_flow;
    port_a.h_outflow=h[1];
    port_b.h_outflow=h[n];

    state_a=Medium.setState_ph(port_a.p,inStream(port_a.h_outflow));
    state_b=Medium.setState_ph(port_b.p,inStream(port_b.h_outflow));

    // Uniform pressure
    p[2:n]=fill(p[1],n-1);

    port_a.p=p[1];
    port_b.p=p[n];

    connect(heatPorts,heatTransfer.heatPorts);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end RefFlow1D_UniformPressure;

  model Compressor "Chiller centrifugal compressor"

    import Modelica.Units.Conversions.to_degC;

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialCompressor;

    Modelica.Blocks.Interfaces.RealInput gamma "Inlet guide vane position" annotation (Placement(
          transformation(extent={{-2,-110},{16,-92}}), iconTransformation(
          extent={{-9,-9},{9,9}},
          rotation=90,
          origin={1,-101})));
    Modelica.Blocks.Interfaces.RealOutput Q_dot_m "Motor heat loss";

    // Map coefficients
    constant Real c[:]={1.4354,0.0054572,-0.0030135,-0.017697,0.000005653};
    constant Real b[:]={7.2058,0.8,0.003};
    //constant Real a[:]={-0.26524,7.1149,-23.415,0.04173,-0.00089576};
    constant Real a[:]={0.0152957946,6.06679198,-17.8814626,0.00936346555,0.000159419745};
    constant SI.Efficiency eta_em=0.75 "Motor electro-mechanical efficiency";

  protected
    constant Real RLA_limit=110;
    SI.MassFlowRate m_flow_max(start=m_flow_init);
    SI.VolumeFlowRate V_dot;
    SI.SpecificEnergy W_p "Ploytropic work";
    SI.Efficiency eta_p(min=0.01,max=1.0) "Polytropic efficiency";
    Real RLA;
    Medium.Density rho_dis=Medium.density(state_dis);
    Medium.Density rho_suc=Medium.density(state_suc);

  equation

    V_dot=m_flow/rho_suc;
    W_p=(p_dis/rho_dis-p_suc/rho_suc)*log(p_dis/p_suc)/log((p_dis/rho_dis)/(p_suc/rho_suc));
    Pwr=m_flow*W_p/(eta_p*eta_em);
    Q_loss=Pwr*(1-eta_em);
    h_dis=h_suc+Pwr*eta_em/m_flow;
    m_flow=gamma*m_flow_max;

    /* Mass balance */
    port_a.m_flow + port_b.m_flow = 0;

    /*   Mass flow rate map  */
    m_flow_max=c[1]+c[2]*p_suc/1e3+c[3]*p_dis/1e3+c[4]*to_degC(T_suc)+c[5]*p_suc/1e3*p_dis/1e3;

    /*  Polytropic efficiency map  */
    eta_p=a[1]+a[2]*V_dot+a[3]*V_dot^2+a[4]*W_p/1e3+a[5]*(W_p/1e3)^2; // Wp in [kJ/kg] in the map

    /*  RLA map  */
    RLA=b[1]+b[2]*Pwr/1e3+b[3]*(Pwr/1e3)^2;

    Q_dot_m=Q_loss;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Compressor;

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

  model TXV "Thermo-static expansion valve"

    extends DynamicVCC.Components.Units.MassFlowDevices.BaseClasses.PartialValve;

    Modelica.Blocks.Interfaces.RealInput p_b "Bulb pressure" annotation (Placement(
          transformation(extent={{-8,-66},{24,-34}}), iconTransformation(
          extent={{-11,-11},{11,11}},
          rotation=90,
          origin={-1,-43})));
    Modelica.Blocks.Interfaces.RealInput Q_dot_m "Motor cooling heat";

    import Modelica.Constants.eps;

  protected
    parameter SI.AbsolutePressure p_min=80e3;
    parameter SI.Area A_max=250e-6;
    parameter Real k_spring=42e-5 "Spring constant";
    parameter Real coef[2]={0.01764338, -0.24770373} "Coefficients for Area with lift";
    parameter Real lift_max=0.0245092036131709;
    SI.Pressure delta_p "Pressure difference bulb and suction refrigerant";
    Medium.MassFlowRate m_dot_cl "chiller motor cooling line";
    Medium.MassFlowRate m_dot_v "valve flow rate";
    SI.Area Avalve;
    Real lift;
    Medium.Density rho_in=Medium.density(state_a);

  equation
    Cd=0.4;
    delta_p=max(eps,(p_b-port_b.p-p_min));
    lift=min(lift_max,k_spring*delta_p/1e3);
    Avalve=max(eps,min(A_max,coef[1]*lift+coef[2]*lift^2));
    m_dot_v=Cd*Avalve*sqrt(rho_in*max(dp,eps));
    m_dot_cl=Cd*8e-6*sqrt(rho_in*max(dp,eps));
    m_flow=m_dot_v+m_dot_cl;
    (port_b.h_outflow-port_a.h_outflow)*m_flow=Q_dot_m;

    port_a.h_outflow=inStream(port_a.h_outflow);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end TXV;

  function vaneAction

     input Real abserr "absolute error";
     output Real actionTime;
  protected
     parameter Real deadband=0.01;
     parameter Real modlimit=1.06;
     parameter Real maxstep=10;
     parameter Real waittime=15;

  algorithm
     if abserr<=deadband then
       actionTime:=0;
     elseif abserr>deadband and abserr<=modlimit then
       actionTime:=(abserr - deadband)/(modlimit - deadband)*maxstep;
     else
       actionTime:=waittime;
     end if;

  end vaneAction;

  function controlSignal

    input Boolean open;
    input Real gamma_pre;
    input Real actionTime_pre;
    output Real gamma;
    output Real actionTime;

  protected
    parameter Real gamma_max=1;
    parameter Real gamma_min=0.05;
    parameter Real upmax=1/150;
    parameter Real dnmax=1/38;//1/45;
    parameter Real waitTime=15;
  algorithm
    if actionTime_pre<1 then
      if open and gamma_pre<gamma_max then
        gamma:=gamma_pre + upmax*actionTime_pre;
      elseif (not open) and gamma_pre>gamma_min then
        gamma:=gamma_pre - dnmax*actionTime_pre;
      else
        gamma:=gamma_pre;
      end if;
      actionTime:= 0;
    else
      if open and gamma_pre<gamma_max then
        gamma:=gamma_pre + upmax;
      elseif (not open) and gamma_pre>gamma_min then
        gamma:=gamma_pre - dnmax;
      else
        gamma:=gamma_pre;
      end if;
      actionTime:=actionTime_pre - 1;
    end if;

  end controlSignal;

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

  model Cycle
    // Calibrating 10 parameters
    Modelica.Blocks.Interfaces.RealOutput y[7] "Outputs of pressures, water exit temperatures, superheat, subcooling, power";
    Modelica.Blocks.Interfaces.RealOutput y_mea[7];

    extends DynamicVCC.Examples.Chiller.HeatExchangers(
      final alpha_f_cond=u[1],
      final alpha_tp_cond=u[2],
      final alpha_g_cond=u[3],
      final alpha_w_cond=u[4],
      final C_metalWall_cond=u[5],
      final alpha_f_evap=1e5,
      final alpha_tp_evap=u[6],
      final alpha_g_evap=u[7],
      final alpha_w_evap=u[8],
      final C_metalWall_evap=u[9]);


    //parameter Real u[10]={113124.022, 111309.067, 4466.849, 124994.011, 108085.516, 12971.468, 882566.944, 32128.593, 72140.224, 142.217};
    parameter Real u[10]={121000.067,114288.143,99973.9913,122046.029,84188.9236,13776.0744,794799.445,17782.3532,73620.7346,144.237819};

    /************** Components **************/

    parameter SI.Power Pwr_init=80e3;
    parameter SI.Temperature Tb_init=285;
    parameter Real gamma_init=0.8624;
    parameter Real C_b=u[10];

    Compressor compressor(
    redeclare package Medium=Medium_1,
    m_flow_init=m_flow_init,
    Pwr_init=Pwr_init,
    h_dis_init=h_init_cond[1],
    h_suc_init=h_init_evap[Ncell],
    p_dis_init=p_init_cond,
    p_suc_init=p_init_evap) annotation (Placement(transformation(
          extent={{31,-31},{-31,31}},
          rotation=-90,
          origin={111,3})));

    Bulb bulb(
    final C=C_b,
    T_init=Tb_init,
    SteadyState_init=false) annotation (Placement(transformation(extent={{-154,-54},{-134,-34}})));

    Modelica.Fluid.Sensors.Temperature temperature_suc(
    redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{90,-96},{66,-72}})));

    Controllernew controller(
    gamma_init=gamma_init) annotation (Placement(transformation(extent={{14,-14},{34,6}})));

    Modelica.Fluid.Sources.MassFlowSource_T sourceCond(nPorts=1,
    redeclare package Medium=Medium_2,
    final use_m_flow_in=true,
    final use_T_in=true) annotation (Placement(transformation(extent={{-58,8},{-38,28}})));

    Modelica.Fluid.Sources.FixedBoundary sinkCond(nPorts=1,
    redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{80,64},{60,84}})));

    Modelica.Fluid.Sources.MassFlowSource_T sourceEvap(
    redeclare package Medium=Medium_2,
    final use_m_flow_in=true,
    final use_T_in=true,
      nPorts=1)          annotation (Placement(transformation(extent={{66,-52},{46,-32}})));

    Modelica.Fluid.Sources.FixedBoundary sinkEvap(nPorts=1,
    redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{-70,-50},{-50,-30}})));

    Modelica.Fluid.Sensors.Temperature temperature_supply(
    redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{-34,-30},{-10,-6}})));

    TXV txv(
    redeclare package Medium=Medium_1,
    dp_nominal=500000,
    m_flow_init=m_flow_init,
    m_flow_nominal=2) annotation (Placement(transformation(
          extent={{21,21},{-21,-21}},
          rotation=90,
          origin={-87,-1})));

    Components.Units.Sensors.T_Superheat subcooling(redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-124,54},{-104,74}})));
    Components.Units.Sensors.T_Superheat superheating(redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{128,-56},{148,-36}})));

    SI.Mass charge;

    Modelica.Fluid.Sensors.Temperature temperature_cond(redeclare package Medium = Medium_2) "Condenser water exit temperature" annotation (Placement(transformation(extent={{132,68},{156,92}})));
    Modelica.Blocks.Sources.CombiTimeTable BC(
      tableOnFile=true,
      tableName="BC",
      fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC.mat",
      columns=2:10,
      smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
      timeScale=1)                                                "Cycle boundary conditions" annotation (Placement(transformation(extent={{-192,60},{-172,80}})));
    Modelica.Blocks.Sources.CombiTimeTable Mea(
      tableOnFile=true,
      tableName="Mea",
      fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",
      columns=2:11,
      smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
      timeScale=1)                                                "Measurement" annotation (Placement(transformation(extent={{-192,30},{-172,50}})));

  initial equation
    //charge=50;

  equation

    charge=condenser.charge+evaporator.charge;


    y[1]=compressor.port_a.p;
    y[2]=compressor.port_b.p;
    y[3]=temperature_supply.T;
    y[4]=temperature_cond.T;
    y[5]=superheating.T;
    y[6]=subcooling.T;
    y[7]=compressor.Pwr;

    y_mea[1]=Mea.y[1];
    y_mea[2]=Mea.y[2];
    y_mea[3]=Mea.y[3];
    y_mea[4]=Mea.y[4];
    y_mea[5]=Mea.y[7];
    y_mea[6]=Mea.y[8];
    y_mea[7]=Mea.y[10];



    controller.u2=BC.y[9];
    sourceCond.m_flow_in=16.9;
    sourceCond.T_in=BC.y[6];
    sourceEvap.m_flow_in=13.6;
    sourceEvap.T_in=BC.y[5];
    /*
  controller.u2=282.55;
  sourceCond.m_flow_in=16.9;
  sourceCond.T_in=302.95;
  sourceEvap.m_flow_in=13.6;
  sourceEvap.T_in=289.05;
*/

    connect(compressor.Q_dot_m,txv.Q_dot_m);

    connect(evaporator.port_b1, compressor.port_a) annotation (Line(points={{22,-63},{110,-63},{110,-30},{111,-30},{111,-28}}, color={0,127,255}));
    connect(temperature_suc.T, bulb.u) annotation (Line(points={{69.6,-84},{36,-84},{36,-92},{-166,-92},{-166,-44},{-156,-44}}, color={0,0,127}));
    connect(sourceCond.ports[1], condenser.port_a2) annotation (Line(points={{-38,18},{-36,18},{-36,49.98},{-21.58,49.98}}, color={0,127,255}));
    connect(condenser.port_b2, sinkCond.ports[1]) annotation (Line(points={{20,75.6},{20,74},{60,74}},           color={0,127,255}));
    connect(sinkEvap.ports[1], evaporator.port_b2) annotation (Line(points={{-50,-40},{-50,-50.4},{-20,-50.4}}, color={0,127,255}));
    connect(temperature_supply.T, controller.u1) annotation (Line(points={{-13.6,-18},{-2,-18},{-2,2},{12,2}}, color={0,0,127}));
    connect(controller.y, compressor.gamma) annotation (Line(points={{35,-4},{74,-4},{74,3.31},{84.9083,3.31}}, color={0,0,127}));
    connect(temperature_supply.port, evaporator.port_b2) annotation (Line(points={{-22,-30},{-40,-30},{-40,-50.4},{-20,-50.4}}, color={0,127,255}));
    connect(bulb.y, txv.p_b) annotation (Line(points={{-133,-44},{-114,-44},{-114,-0.79},{-96.03,-0.79}}, color={0,0,127}));
    connect(subcooling.port, txv.port_a) annotation (Line(points={{-114,54},{-114,26},{-87,26},{-87,20}},        color={0,127,255}));
    connect(superheating.port, compressor.port_a) annotation (Line(points={{138,-56},{110,-56},{110,-30},{111,-30},{111,-28}}, color={0,127,255}));
    connect(condenser.port_b2, temperature_cond.port) annotation (Line(points={{20,75.6},{20,74},{54,74},{54,90},{126,90},{126,62},{144,62},{144,68}},
                                                                                                                             color={0,127,255}));

    connect(condenser.port_b1, txv.port_a) annotation (Line(points={{-22,63},{-22,62},{-87,62},{-87,20}}, color={0,127,255}));
    connect(condenser.port_a1, compressor.port_b) annotation (Line(points={{20,63},{20,62},{54,62},{54,40},{111,40},{111,34}}, color={0,127,255}));
    connect(txv.port_b, evaporator.port_a1) annotation (Line(points={{-87,-22},{-87,-63},{-20,-63}}, color={0,127,255}));
    connect(evaporator.port_b1, temperature_suc.port) annotation (Line(points={{22,-63},{98,-63},{98,-102},{78,-102},{78,-96}}, color={0,127,255}));
    connect(sourceEvap.ports[1], evaporator.port_a2) annotation (Line(points={{46,-42},{40,-42},{40,-76.02},{21.58,-76.02}}, color={0,127,255}));
    annotation (Diagram(coordinateSystem(extent={{-200,-100},{180,100}})), Icon(coordinateSystem(extent={{-200,-100},{180,100}}), graphics={
          Rectangle(
            extent={{-34,42},{36,76}},
            lineColor={0,0,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-36,-96},{34,-62}},
            lineColor={0,0,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{148,20},{168,20},{188,-20},{128,-20},{148,20}},
            lineColor={0,0,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-160,4},{-170,24},{-150,24},{-150,24},{-160,4}},
            lineColor={0,0,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-160,4},{-170,-16},{-150,-16},{-150,-16},{-160,4}},
            lineColor={0,0,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Line(
            points={{36,60},{160,60},{160,20}},
            color={0,0,0},
            thickness=1),
          Line(
            points={{34,-80},{162,-80},{162,-20}},
            color={0,0,0},
            thickness=1),
          Line(
            points={{-4,31},{32,31},{32,-95}},
            color={0,0,0},
            thickness=1,
            origin={-129,28},
            rotation=90),
          Line(
            points={{-91,44},{33,44},{33,-20}},
            color={0,0,0},
            thickness=1,
            origin={-127,-36},
            rotation=180)}),
      experiment(
        StartTime=2000,
        StopTime=39500,
        Interval=10,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentFlags(
        Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
        Evaluate=false,
        OutputCPUtime=true,
        OutputFlatModelica=false));
  end Cycle;

  model Controllernew
   extends Modelica.Blocks.Interfaces.SI2SO; //u1=Twater, u2=Tset, y=gamma

   parameter Real starttime=0;
   parameter Real gamma_init;

  protected
    parameter Real gamma_max=1;
    parameter Real gamma_min=0.05;
    Real actionTime;
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
      actionTime:=vaneAction(abs(error));
    end when;

    when sample(starttime+0.01,sampletime_vane) then
        (gamma,actionTime):=controlSignal(open,pre(gamma),pre(actionTime));
      end when;

    y:=max(min(gamma_max,gamma),gamma_min);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end Controllernew;

  package Test "Test component models of chiller"
    extends Modelica.Icons.ExamplesPackage;

    model TestCompressor
      extends Modelica.Icons.Example;

      inner DynamicVCC.Components.System system(
      p_max=12e5,
      p_min=2e5,
      h_max=5e5,
      h_min=1e5,
      T_max=340,
      T_min=260,
      EnableReverseFlow=false);

      replaceable package Medium=Modelica.Media.R134a.R134a_ph;

      DynamicVCC.Examples.Chiller.Compressor compressor(redeclare package Medium = Medium);

      // Suction
      Modelica.Fluid.Sources.Boundary_ph source(
      redeclare package Medium=Medium,
      nPorts=1,
      use_p_in=true,
      use_h_in=true);

      // Discharge
      Modelica.Fluid.Sources.Boundary_ph sink(
      redeclare package Medium=Medium,
      nPorts=1,
      use_p_in=true);

      Modelica.Blocks.Sources.CombiTimeTable BC_Cond(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Cond",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Cond.mat",columns=2:5);
      Modelica.Blocks.Sources.CombiTimeTable BC_Evap(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Evap",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Evap.mat",columns=2:5);
      Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:9);

    equation

      compressor.gamma=1;

      source.p_in=Mea.y[1];
      source.h_in=Mea.y[5];
      sink.p_in=Mea.y[2];

      connect(source.ports[1],compressor.port_a);
      connect(compressor.port_b,sink.ports[1]);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StartTime=2000,
          StopTime=5000,
          __Dymola_Algorithm="Dassl"));
    end TestCompressor;

    model TestController
      extends Modelica.Icons.Example;

      DynamicVCC.Examples.Chiller.Controllernew controller(gamma_init=0.9096);

      Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:10);
      Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC.mat",columns=2:10);

    equation
      controller.u1=Mea.y[3];
      controller.u2=BC.y[9];

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StartTime=2000,
          StopTime=50000,
          __Dymola_Algorithm="Dassl"));
    end TestController;

    model TestTXV
      extends Modelica.Icons.Example;

      inner DynamicVCC.Components.System system(
        p_max=12e5,
        p_min=2e5,
        h_max=5e5,
        h_min=1e5,
        T_max=340,
        T_min=260,
        EnableReverseFlow=false);

      replaceable package Medium=Modelica.Media.R134a.R134a_ph;

      DynamicVCC.Examples.Chiller.Bulb bulb(C=50);

      DynamicVCC.Examples.Chiller.TXV txv(
        redeclare package Medium = Medium,
        dp_nominal=5e5,
        m_flow_init=2,
        m_flow_nominal=2);

      Medium.ThermodynamicState state_suc; //suction state

      //Source
      Modelica.Fluid.Sources.Boundary_ph source(
      redeclare package Medium=Medium,
      nPorts=1,
      use_p_in=true,
      use_h_in=true);

      //Sink
      Modelica.Fluid.Sources.Boundary_ph sink(
      redeclare package Medium=Medium,
      nPorts=1,
      use_p_in=true);

      Modelica.Blocks.Sources.CombiTimeTable BC_Cond(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Cond",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Cond.mat",columns=2:5);
      Modelica.Blocks.Sources.CombiTimeTable BC_Evap(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC_Evap",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC_Evap.mat",columns=2:5);
      Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:9);

    equation
      state_suc=Medium.setState_ph(Mea.y[1],Mea.y[5]);
      bulb.u=Medium.temperature(state_suc);
      txv.Q_dot_m=0;

      source.p_in=Mea.y[2];
      source.h_in=Mea.y[6];
      sink.p_in=Mea.y[1];

      connect(bulb.y,txv.p_b);
      connect(source.ports[1],txv.port_a);
      connect(txv.port_b,sink.ports[1]);

      annotation (experiment(
          StartTime=2000,
          StopTime=5000,
          Tolerance=0.001,
          __Dymola_Algorithm="Dassl"));
    end TestTXV;

    model Cycle_BO "Cycle model for Bayesian Optimization"

      Modelica.Blocks.Interfaces.RealOutput y[6] "Outputs of pressures, water exit temperatures, superheat, subcooling";

      parameter Real u[10]={250423.409000000,870705.033000000,339941.861000000,155099.166000000,1464890.11000000,272490.187000000,1050636.90000000,216781.467000000,363297.551000000,98.0396584000000};

      //{5.97070338e+05, 2.86528623e+06, 2.34213860e+06, 1.15112818e+05,6.44915651e+05, 1.15699143e+06, 6.49488319e+05, 7.24785181e+04,2.57169763e+05, 7.01628514e+01};

      // Note: heat transfer coefficients are actually UA, area of HX model is set to 1.
      // And metal mass is actually heat capacity, specific heat capacity of metal material is set to 1.
      parameter SI.CoefficientOfHeatTransfer alpha_cond=u[1];
      parameter SI.CoefficientOfHeatTransfer alpha_f_c=u[2];
      parameter SI.CoefficientOfHeatTransfer alpha_g_c=u[3];
      parameter SI.CoefficientOfHeatTransfer alpha_evap=u[4];
      parameter SI.CoefficientOfHeatTransfer alpha_g_e=u[5];
      parameter SI.CoefficientOfHeatTransfer alpha_w_c=u[6];
      parameter SI.CoefficientOfHeatTransfer alpha_w_e=u[7];
      parameter SI.Mass Mw_c=u[8];
      parameter SI.Mass Mw_e=u[9];
      parameter Real C_b=u[10] "Bulb time constant";

      inner Components.System system(
      p_max=12e5,
      p_min=2e5,
      h_max=5e5,
      h_min=1e5,
      T_max=340,
      T_min=260,
      EnableReverseFlow=EnableReverseFlow) annotation (Placement(transformation(extent={{-162,48},{-126,84}})));

      package Medium_CP=DynamicVCC.Media.CoolProp.R134a;
      //package Medium_NN=DynamicVCC.Media.R134a_NN;

      package Medium_1=Medium_CP;

      package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater (
      cp_const=4186.8,
      d_const=995,
      lambda_const=0.6);

      /************** Initial conditions **************/
      parameter SI.Pressure p_init_c[Ncell]=fill(9.0556e+05,Ncell);
      parameter SI.SpecificEnthalpy h_init_c[Ncell]={425665,422357.093750000,420374.250000000,419181.812500000,418463.156250000,418029.500000000,417767.562500000,414636.218750000,409043.187500000,399053.312500000,381210.031250000,349339.531250000,292414.468750000,250096.875000000,247136.843750000};
      parameter SI.ThermodynamicTemperature Tt_init_c[Ncell]={309.448150634766,309.167724609375,309.000793457031,308.900909423828,308.840911865234,308.804748535156,308.782958984375,308.891448974609,308.885467529297,308.874816894531,308.855804443359,308.821807861328,308.761108398438,306.224029541016,303.171966552734};
      parameter SI.ThermodynamicTemperature Te_init_c[Ncell]={303.050842285156,304.492584228516,306.432006835938,307.517822265625,308.125732421875,308.466064453125,308.656616210938,308.763305664063,308.772247314453,308.787017822266,308.811492919922,308.852111816406,308.919677734375,309.032379150391,309.221282958984};

      //parameter SI.Pressure p_init_e[Ncell]={381953.343750000,381953.343750000,381953.343750000,381953.343750000,381953.343750000,381953.343750000,381953.343750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000,381297.093750000};
      parameter SI.Pressure p_init_e[Ncell]=fill(381953.34375,Ncell);
      parameter SI.SpecificEnthalpy h_init_e[Ncell]={248823.453125000,255712.937500000,263863.593750000,273506.250000000,284914.093750000,298410.187500000,314376.812500000,333487.500000000,356096.468750000,382844.187500000,402964.093750000,408547.500000000,410089.375000000,410514.781250000,410632.125000000};
      parameter SI.ThermodynamicTemperature Tt_init_e[Ncell]={280.904602050781,280.938842773438,280.979370117188,281.027343750000,281.084045410156,281.151153564453,281.230560302734,281.281097412109,281.393524169922,281.526519775391,283.943664550781,287.634765625000,288.660552978516,288.943908691406,289.022094726563};
      parameter SI.ThermodynamicTemperature Te_init_e[Ncell]={289.045104980469,289.027313232422,288.962921142578,288.729644775391,287.889068603516,286.771575927734,285.826995849609,285.028564453125,284.361511230469,283.797668457031,283.321044921875,282.918182373047,282.577667236328,282.289825439453,282.046539306641};

      parameter SI.MassFlowRate m_flows_init[Ncell+1]=fill(2.396393,Ncell+1);
      parameter Boolean SteadyState_init=false;
      parameter SI.Pressure dp_init=5e5;
      parameter SI.Power Pwr_init=80e3;
      parameter Real gamma_init=0.9;
      parameter Real Tb_init=284;

      /************** Numerical **************/
      parameter Integer Ncell=15;

      import DynamicVCC.Components.Types.ModelStructure;
      parameter Boolean EnableReverseFlow=true;
      parameter Boolean useLumpedPressure=true;
      parameter Boolean use_I_flows=true;
      parameter ModelStructure modelStructure=ModelStructure.av_vb;
      import DynamicVCC.Components.Types.DifferentialState;
      parameter DifferentialState differentialState=DifferentialState.pdh;

       /************** Heat transfer **************/
      replaceable model HeatTransfer_1_c=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
      final alpha_f=alpha_f_c,
      final alpha_tp=alpha_cond,
      final alpha_g=alpha_g_c);

    /*
  replaceable model LiquidPhase_e =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=3e3,
  final alpha_cst=alpha_f_e);

  replaceable model VaporPhase_e =DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Constant (
  final alpha0=3e3,
  final alpha_cst=alpha_g_e);


  replaceable model HeatTransfer_1_e=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.HeatTransferPhaseZones (
  redeclare model LiquidZone=LiquidPhase_e,
  redeclare model VaporZone=VaporPhase_e,
  redeclare model TwoPhaseZone=Evaporation);
  */

      replaceable model HeatTransfer_1_e=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
      final alpha_f=1e4,
      final alpha_tp=alpha_evap,
      final alpha_g=alpha_g_e);

      replaceable model HeatTransfer_2_c = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.ConstHTC (
      final alpha0=alpha_w_c);

      replaceable model HeatTransfer_2_e = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.ConstHTC (
      final alpha0=alpha_w_e);
    /*
  replaceable model HeatTransfer_2_c=DynamicVCC.Test.Chiller.HeatTransfer.WaterCond (
  final C_sf=8);

  replaceable model HeatTransfer_2_e=DynamicVCC.Test.Chiller.HeatTransfer.WaterEvap (
  final C_sf=3);
  */

      /************** Friction **************/
      replaceable model Friction_1_c = DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant(f0=0.01);

      replaceable model Friction_1_e = DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant(f0=0.01);

      /************** Components **************/

      replaceable model RefFlow=DynamicVCC.Examples.Chiller.RefFlow1D_UniformPressure;

      Components.Units.HX.ShellTubeHX condenser(
      redeclare final package Medium_1=Medium_1,
      redeclare final package Medium_2=Medium_2,
      redeclare final model RefFlow1D=RefFlow,
      redeclare final model HeatTransfer_1=HeatTransfer_1_c,
      redeclare final model HeatTransfer_2=HeatTransfer_2_c,
      redeclare final model Friction_1=Friction_1_c,
      final Ncell=Ncell,
      As_1=1,
      Ac_1=0.0652,
      L_1=2.4384,
      diameter_1=0.01905,
      As_2=1,
      Ac_2=0.0311,
      L_2=2.4384,
      diameter_2=0.01554,
      M_metalWall=Mw_c,
      cp_metalWall=1,
      modelStructure=modelStructure,
      differentialState=differentialState,
      EnableReverseFlow=EnableReverseFlow,
      useLumpedPressure=useLumpedPressure,
      use_I_flows=use_I_flows,
      SteadyState_init=SteadyState_init,
      p_init=p_init_c,
      h_init=h_init_c,
      Tt_init=Tt_init_c,
      Te_init=Te_init_c,
      m_flows_init=m_flows_init) annotation (Placement(transformation(extent={{28,30},{-28,86}})));

      Components.Units.HX.ShellTubeHX evaporator(
      redeclare final package Medium_1=Medium_1,
      redeclare final package Medium_2=Medium_2,
      redeclare final model RefFlow1D=RefFlow,
      redeclare final model HeatTransfer_1=HeatTransfer_1_e,
      redeclare final model HeatTransfer_2=HeatTransfer_2_e,
      redeclare final model Friction_1=Friction_1_e,
      final Ncell=Ncell,
      As_1=1,
      Ac_1=0.07627,
      L_1=2.4384,
      diameter_1=0.0196,
      As_2=1,
      Ac_2=0.03,
      L_2=2.4384,
      diameter_2=0.01606,
      M_metalWall=Mw_e,
      cp_metalWall=1,
      modelStructure=modelStructure,
      differentialState=differentialState,
      EnableReverseFlow=EnableReverseFlow,
      useLumpedPressure=useLumpedPressure,
      use_I_flows=use_I_flows,
      SteadyState_init=SteadyState_init,
      p_init=p_init_e,
      h_init=h_init_e,
      Tt_init=Tt_init_e,
      Te_init=Te_init_e,
      m_flows_init=m_flows_init) annotation (Placement(transformation(extent={{-28,-86},{28,-30}})));

      Compressor compressor(
      redeclare package Medium=Medium_1,
      m_flow_init=m_flows_init[1],
      Pwr_init=Pwr_init) annotation (Placement(transformation(
            extent={{31,-31},{-31,31}},
            rotation=-90,
            origin={111,5})));

      Bulb bulb(
      final C=C_b,
      T_init=Tb_init,
      SteadyState_init=SteadyState_init) annotation (Placement(transformation(extent={{-154,-54},{-134,-34}})));

      Modelica.Fluid.Sensors.Temperature temperature_suc(
      redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{86,-92},{62,-68}})));

      Controllernew controller(
      gamma_init=gamma_init) annotation (Placement(transformation(extent={{14,-14},{34,6}})));

      Modelica.Fluid.Sources.MassFlowSource_T sourceCond(nPorts=1,
      redeclare package Medium=Medium_2,
      final use_m_flow_in=true,
      final use_T_in=true) annotation (Placement(transformation(extent={{-58,8},{-38,28}})));

      Modelica.Fluid.Sources.FixedBoundary sinkCond(nPorts=1,
      redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{68,68},{48,88}})));

      Modelica.Fluid.Sources.MassFlowSource_T sourceEvap(nPorts=1,
      redeclare package Medium=Medium_2,
      final use_m_flow_in=true,
      final use_T_in=true) annotation (Placement(transformation(extent={{66,-52},{46,-32}})));

      Modelica.Fluid.Sources.FixedBoundary sinkEvap(nPorts=1,
      redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{-70,-50},{-50,-30}})));

      Modelica.Fluid.Sensors.Temperature temperature_supply(
      redeclare package Medium=Medium_2) annotation (Placement(transformation(extent={{-34,-30},{-10,-6}})));

      Modelica.Blocks.Sources.CombiTimeTable Mea(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="Mea",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/Mea.mat",columns=2:11);
      Modelica.Blocks.Sources.CombiTimeTable BC(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="BC",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/BC.mat",columns=2:10);
      //Modelica.Blocks.Sources.CombiTimeTable Caldata(tableOnFile=true,smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,tableName="data_calibration",fileName="C:/Jiacheng Ma/BoundaryCondition/Chiller/ChillerSimulationData.mat",columns=2:7);

      TXV txv(
      redeclare package Medium=Medium_1,
        dp_nominal=500000,
      m_flow_init=m_flows_init[1],
      m_flow_nominal=2) annotation (Placement(transformation(
            extent={{21,21},{-21,-21}},
            rotation=90,
            origin={-87,-1})));

      Components.Units.Sensors.T_Superheat subcooling(redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{-92,64},{-72,84}})));
      Components.Units.Sensors.T_Superheat superheating(redeclare package Medium=Medium_1) annotation (Placement(transformation(extent={{128,-56},{148,-36}})));

      SI.Mass charge;

      Modelica.Fluid.Sensors.Temperature temperature_cond(redeclare package Medium = Medium_2) "Condenser water exit temperature" annotation (Placement(transformation(extent={{72,70},{96,94}})));
    equation

      charge=condenser.charge+evaporator.charge;
      /*
  alpha_cond=2e4;
  alpha_f_c=1e4;
  alpha_g_c=1e4;
  alpha_w_c=3e4;

  alpha_evap=1e4;
  alpha_f_e=5e3;
  alpha_g_e=5e3;
  alpha_w_e=3e4;
*/

      y[1]=compressor.port_a.p;
      y[2]=compressor.port_b.p;
      y[3]=temperature_supply.T;
      y[4]=temperature_cond.T;
      y[5]=superheating.T;
      y[6]=subcooling.T;

      controller.u2=BC.y[9];
      sourceCond.m_flow_in=16.9;
      sourceCond.T_in=BC.y[6];
      sourceEvap.m_flow_in=13.6;
      sourceEvap.T_in=BC.y[5];
      /*
  controller.u2=282.55;
  sourceCond.m_flow_in=16.9;
  sourceCond.T_in=302.95;
  sourceEvap.m_flow_in=13.6;
  sourceEvap.T_in=289.05;
  */

      connect(compressor.Q_dot_m,txv.Q_dot_m);

      connect(compressor.port_b, condenser.port_a1) annotation (Line(points={{111,36},{110,36},{110,58},{28,58}}, color={0,127,255}));
      connect(evaporator.port_b1, compressor.port_a) annotation (Line(points={{28,-58},{110,-58},{110,-30},{111,-30},{111,-26}}, color={0,127,255}));
      connect(evaporator.port_b1, temperature_suc.port) annotation (Line(points={{28,-58},{94,-58},{94,-98},{74,-98},{74,-92}}, color={0,127,255}));
      connect(temperature_suc.T, bulb.u) annotation (Line(points={{65.6,-80},{36,-80},{36,-92},{-166,-92},{-166,-44},{-156,-44}}, color={0,0,127}));
      connect(sourceCond.ports[1], condenser.port_a2) annotation (Line(points={{-38,18},{-36,18},{-36,40.64},{-27.44,40.64}}, color={0,127,255}));
      connect(condenser.port_b2, sinkCond.ports[1]) annotation (Line(points={{28,74.8},{38,74.8},{38,78},{48,78}}, color={0,127,255}));
      connect(sourceEvap.ports[1], evaporator.port_a2) annotation (Line(points={{46,-42},{36,-42},{36,-75.36},{27.44,-75.36}}, color={0,127,255}));
      connect(sinkEvap.ports[1], evaporator.port_b2) annotation (Line(points={{-50,-40},{-50,-41.2},{-28,-41.2}}, color={0,127,255}));
      connect(temperature_supply.T, controller.u1) annotation (Line(points={{-13.6,-18},{-2,-18},{-2,2},{12,2}}, color={0,0,127}));
      connect(controller.y, compressor.gamma) annotation (Line(points={{35,-4},{74,-4},{74,5.31},{84.9083,5.31}}, color={0,0,127}));
      connect(temperature_supply.port, evaporator.port_b2) annotation (Line(points={{-22,-30},{-40,-30},{-40,-41.2},{-28,-41.2}}, color={0,127,255}));
      connect(condenser.port_b1, txv.port_a) annotation (Line(points={{-28,58},{-87,58},{-87,20}}, color={0,127,255}));
      connect(txv.port_b, evaporator.port_a1) annotation (Line(points={{-87,-22},{-87,-58},{-28,-58}}, color={0,127,255}));
      connect(bulb.y, txv.p_b) annotation (Line(points={{-133,-44},{-114,-44},{-114,-0.79},{-96.03,-0.79}}, color={0,0,127}));
      connect(subcooling.port, txv.port_a) annotation (Line(points={{-82,64},{-68,64},{-68,58},{-87,58},{-87,20}}, color={0,127,255}));
      connect(superheating.port, compressor.port_a) annotation (Line(points={{138,-56},{110,-56},{110,-30},{111,-30},{111,-26}}, color={0,127,255}));
      connect(condenser.port_b2, temperature_cond.port) annotation (Line(points={{28,74.8},{40,74.8},{40,62},{84,62},{84,70}}, color={0,127,255}));
                         annotation (Placement(transformation(
            extent={{-25,-25},{25,25}},
            rotation=-90,
            origin={-87,1})),
                  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},{180,100}}), graphics={
            Rectangle(
              extent={{-40,44},{30,78}},
              lineColor={0,0,0},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-36,-96},{34,-62}},
              lineColor={0,0,0},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{82,20},{102,20},{122,-20},{62,-20},{82,20}},
              lineColor={0,0,0},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-90,4},{-100,24},{-80,24},{-80,24},{-90,4}},
              lineColor={0,0,0},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-90,4},{-100,-16},{-80,-16},{-80,-16},{-90,4}},
              lineColor={0,0,0},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid),
            Line(
              points={{32,60},{94,60},{94,20}},
              color={0,0,0},
              thickness=1),
            Line(
              points={{34,-80},{96,-80},{96,-20}},
              color={0,0,0},
              thickness=1),
            Line(
              points={{-4,31},{32,31},{32,-19}},
              color={0,0,0},
              thickness=1,
              origin={-59,28},
              rotation=90),
            Line(
              points={{-19,44},{33,44},{33,-20}},
              color={0,0,0},
              thickness=1,
              origin={-57,-36},
              rotation=180)}),                                                                       Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},{180,100}})),
        experiment(
          StartTime=2000,
          StopTime=50000,
          Interval=10,
          Tolerance=0.001,
          __Dymola_Algorithm="Cvode"));
    end Cycle_BO;
  end Test;

  package HeatTransfer "Heat transfer correlations for chiller"
    extends Modelica.Icons.VariantsPackage;
    model SinglePhase
      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.PartialHeatTransferCorrelation;

      parameter Real C_sf=1 "Surface enhancement factor";

    equation

      Nu=C_sf*sqrt(Re).*Pr;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SinglePhase;

    model WaterEvap "Evaporator water-side heat transfer correlation"
      extends DynamicVCC.Examples.Chiller.HeatTransfer.IncompressibleCorrelation;

      parameter Real C_sf=1 "Surface enhancement factor";

    protected
      SI.CoefficientOfFriction f[n];

    equation
      f={-0.5087+0.2768*log10(Re[i])-0.0339*log10(Re[i])*log10(Re[i]) for i in 1:n};
      Nu={C_sf*(f[i]/8*(Re[i]-1000)*Pr[i])/(1.07+12.7*sqrt(f[i]/8)*(Pr[i]^0.67-1)) for i in 1:n};
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end WaterEvap;

    model WaterCond "Condenser water-side heat transfer correlation"
      extends DynamicVCC.Examples.Chiller.HeatTransfer.IncompressibleCorrelation;

      parameter Real C_sf=1 "Surface enhancement factor";

    protected
      SI.CoefficientOfFriction f[n];

    equation
      f={1/(0.79*log(Re[i])-1.64)^2 for i in 1:n};
      Nu={C_sf*f[i]*Re[i]*Pr[i]/8/(1.07+12.7*sqrt(f[i]/8)*(Pr[i]^0.67-1)) for i in 1:n};
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end WaterCond;

    partial model IncompressibleCorrelation
      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.PartialFlowHeatTransfer;

       // Fluid properties
      Medium.Density rho[n]=Medium.density(states);
      Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states);
      Medium.PrandtlNumber Pr[n]=Medium.prandtlNumber(states);
      Medium.ThermalConductivity k[n]=Medium.thermalConductivity(states);

      // Variables
      SI.ReynoldsNumber Re[n];
      SI.NusseltNumber Nu[n];

    equation
      Re=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(v,rho,mu,dimension);
      Nu=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.NusseltNumber(alphas,dimension,k);
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end IncompressibleCorrelation;
  end HeatTransfer;

  partial model HeatExchangers "Condenser and evapoator models"
    inner DynamicVCC.Components.System system(
    redeclare package Medium=Medium_1,
    T_max=330,
    T_min=250,
    m_flow_init=m_flow_init,
    m_flow_nominal=2.396,
    massDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
    energyDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
    momentumDynamics=DynamicVCC.Components.Types.Dynamics.Fixed_init,
    enableReverseFlow=true);

    package Medium_1=DynamicVCC.Media.R134a_NN;
    package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater (
    cp_const=4186.8,
    d_const=995);

     /******** Initial conditions *********/

    // condenser
    parameter SI.Pressure p_init_cond=9.1e5;
    parameter SI.SpecificEnthalpy h_init_cond[Ncell]={419460.687500000,417711.625000000,415584.750000000,412807.937500000,409182.625000000,404449.531250000,398270.125000000,390267.375000000,379819.187500000,366178.343750000,348369.187500000,325118.062500000,294762.031250000,255130.015625000,243951.390625000};
    parameter SI.ThermodynamicTemperature Tt_init_cond[Ncell]={310.370910644531,309.052246093750,309.027801513672,309.015594482422,308.999694824219,308.978881835938,308.951751708984,308.910339355469,308.864440917969,308.804534912109,308.726318359375,308.624206542969,308.490875244141,308.316802978516,304.357025146484};
    parameter SI.ThermodynamicTemperature Te_init_cond[Ncell]={303.328430175781,304.670043945313,305.697631835938,306.484741210938,307.087585449219,307.549377441406,307.903045654297,308.173950195313,308.383148193359,308.543365478516,308.666107177734,308.760101318359,308.832092285156,308.891296386719,309.289245605469};

    parameter SI.Pressure p_init_evap=4e5;
    parameter SI.SpecificEnthalpy h_init_evap[Ncell]={252861.890625000,253763.671875000,254986.734375000,256645.515625000,258895.265625000,261946.515625000,266084.812500000,271703.937500000,279324.937500000,289661,303679.468750000,322692.187500000,348478.437500000,383451.437500000,406479.062500000};
    parameter SI.ThermodynamicTemperature Tt_init_evap[Ncell]={282.494842529297,282.499603271484,282.506042480469,282.514801025391,282.526641845703,282.542755126953,282.564575195313,282.593536376953,282.633728027344,282.688232421875,282.762145996094,282.862396240234,282.998382568359,283.182800292969,285.824615478516};
    parameter SI.ThermodynamicTemperature Te_init_evap[Ncell]={288.081329345703,286.610168457031,285.525451660156,284.725646972656,284.135955810547,283.701171875000,283.380584716797,283.144195556641,282.970123291016,282.841766357422,282.747131347656,282.677337646484,282.625915527344,282.587982177734,282.559997558594};

    parameter SI.MassFlowRate m_flow_init=2.396;

    /********** Numerics **********/
    parameter Integer Ncell=15;
    import DynamicVCC.Components.Types.ModelStructure;
    import DynamicVCC.Components.Types.DifferentialState;
    parameter ModelStructure modelStructure=ModelStructure.av_vb;
    parameter DifferentialState differentialState=DifferentialState.pdh;
    parameter Boolean useLumpedPressure=true;

    parameter SI.CoefficientOfHeatTransfer alpha_f_cond;
    parameter SI.CoefficientOfHeatTransfer alpha_tp_cond;
    parameter SI.CoefficientOfHeatTransfer alpha_g_cond;
    parameter SI.CoefficientOfHeatTransfer alpha_w_cond;
    parameter SI.HeatCapacity C_metalWall_cond;
    parameter SI.CoefficientOfHeatTransfer alpha_f_evap;
    parameter SI.CoefficientOfHeatTransfer alpha_tp_evap;
    parameter SI.CoefficientOfHeatTransfer alpha_g_evap;
    parameter SI.CoefficientOfHeatTransfer alpha_w_evap;
    parameter SI.HeatCapacity C_metalWall_evap;

     /******** Heat transfer ********/
    replaceable model HeatTransfer_1_cond = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
    final alpha_f=alpha_f_cond,
    final alpha_tp=alpha_tp_cond,
    final alpha_g=alpha_g_cond);

    replaceable model HeatTransfer_2_cond = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    alpha0=alpha_w_cond);

    replaceable model HeatTransfer_1_evap = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowPhaseChange (
    final alpha_f=alpha_f_evap,
    final alpha_tp=alpha_tp_evap,
    final alpha_g=alpha_g_evap);

    replaceable model HeatTransfer_2_evap = DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer (
    alpha0=alpha_w_evap);

    /******** Friction ********/
    replaceable model FlowModel_1_cond = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
      final lambda0=0.01);

    replaceable model FlowModel_1_evap = DynamicVCC.Components.Pipes.BaseClasses.FlowModels.ConstantFrictionFlow (
      final lambda0=0.01);

    /******** Slip Ratio ********/
    replaceable model SlipRatio=DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.Homogeneous;

    /******** Heat exchangers ********/

    Components.Units.HX.ShellTubeHX condenser(
    redeclare final package Medium_1=Medium_1,
    redeclare final package Medium_2=Medium_2,
    redeclare final model HeatTransfer_1=HeatTransfer_1_cond,
    redeclare final model HeatTransfer_2=HeatTransfer_2_cond,
    redeclare final model SlipRatio=SlipRatio,
    redeclare final model FlowModel_1=FlowModel_1_cond,
    final Ncell=Ncell,
    final C_tube=C_metalWall_cond,
    As_1=23.9328,
    Ac_1=0.0652,
    L_1=2.4384,
    diameter_1=0.01905,
    As_2=19.5231,
    Ac_2=0.0311,
    L_2=2.4384,
    diameter_2=0.01554,
    M_metalWall=340.6385,
    cp_metalWall=385,
    modelStructure=modelStructure,
    differentialState=differentialState,
    useLumpedPressure=useLumpedPressure,
    p_a_start=p_init_cond,
    h_init=h_init_cond,
    Tt_init=Tt_init_cond,
    Te_init=Te_init_cond) annotation (Placement(transformation(extent={{20,42},{-22,84}})));


    Components.Units.HX.ShellTubeHX evaporator(
    redeclare final package Medium_1=Medium_1,
    redeclare final package Medium_2=Medium_2,
    redeclare final model HeatTransfer_1=HeatTransfer_1_evap,
    redeclare final model HeatTransfer_2=HeatTransfer_2_evap,
    redeclare final model SlipRatio=SlipRatio,
    redeclare final model FlowModel_1=FlowModel_1_evap,
    final Ncell=Ncell,
    final C_tube=C_metalWall_evap,
    As_1=22.3716,
    Ac_1=0.07627,
    L_1=2.4384,
    diameter_1=0.0196,
    As_2=18.331,
    Ac_2=0.03,
    L_2=2.4384,
    diameter_2=0.01606,
    M_metalWall=321.782954,
    cp_metalWall=385,
    modelStructure=modelStructure,
    differentialState=differentialState,
    useLumpedPressure=useLumpedPressure,
    p_a_start=p_init_evap,
    h_init=h_init_evap,
    Tt_init=Tt_init_evap,
    Te_init=Te_init_evap) annotation (Placement(transformation(extent={{-20,-84},{22,-42}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{140,100}})),
                                                                   Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-100},{140,100}})));
  end HeatExchangers;
end Chiller;
