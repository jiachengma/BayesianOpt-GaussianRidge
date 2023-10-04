within DynamicVCC.Components.Units.HX;
package FrostDefrostHX "Heat exchanger model switching between frosting and defrosting operations"
  extends Modelica.Icons.VariantsPackage;
  model FinTubeHX_FrostDefrost "Switching between frost formation and melting"

    input Real frostmode;

    input SI.Temperature T_amb;

    // extending refrigerant flow models and tube wall models
    extends DynamicVCC.Components.Units.HX.BaseClasses.PartialHX(final C_metalWall=C_FinTube);

    replaceable package Medium_2=Modelica.Media.Air.MoistAir;

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
        DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_ConstCoefficient
    "Air side heat transfer model";

    replaceable model FreeConvection =
      DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer_old.AirCoilHeatTransfer_FreeConvection
    "Air side free convection heat transfer model";

    replaceable model Friction_2 =
        DynamicVCC.Components.Pipes.BaseClasses.Friction.AirCoilDP_ConstFactor
    "Air side frictional pressure drop model";

    // Air flow
    replaceable model AirFlow=DynamicVCC.Components.Pipes.MoistAirCrossFlow;

    AirFlow airFlow(redeclare final package Medium=Medium_2,
    final Ncell=Ncell,
    redeclare final model HeatTransfer=HeatTransfer_2,
    redeclare final model FreeConvection=FreeConvection,
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
    DynamicVCC.Interfaces.FluidPorts_a ports_a2[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
          extent={{-8,-32},{8,32}},
          rotation=90,
          origin={0,-68})));

    DynamicVCC.Interfaces.FluidPorts_b ports_b2[Ncell](redeclare each package Medium = Medium_2) annotation (Placement(transformation(extent={{48,-14},{68,6}}), iconTransformation(
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
      each final T_amb=T_amb);

    replaceable model FrostFormation =
        DynamicVCC.Components.Units.HX.FrostDefrostHX.FrostGrowth_Lee;

    FrostFormation frostFormation[Ncell](
    each As=As_2/Ncell,
    final x_f=x_f,
    final rho_f=rho_f,
    final m_flow_dehumid=airFlow.m_flows_dehumid,
    final T_t=metalWall.heatPorts_a.T);

    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatSink_metalWall[Ncell];

    SI.Mass M_f "Frost mass";

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

    M_f=sum(As_2/Ncell*x_f.*rho_f);

    heatSink_metalWall.Q_flow=Modelica.Fluid.Utilities.regStep(frostmode-0.5,-frostFormation.Q_flow,-frostMelt.Q_flow,1e-3);
    //-frostmode*frostFormation.Q_flow-(1.0-frostmode)*frostMelt.Q_flow;

    connect(ports_a2,airFlow.ports_a);
    connect(ports_b2,airFlow.ports_b);
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

  model FrostGrowth_Lee "(Lee, 1997)"

    import DynamicVCC.Components.Units.HX.FrostDefrostHX.Utilities.waterVaporDensity_sat;
    import thermalConductivity=DynamicVCC.Components.Units.HX.FrostDefrostHX.Utilities.thermalConductivity_Lee;
    import Modelica.Constants.T_zero;
    import Modelica.Fluid.Utilities.regStep;
    import Modelica.Math.acosh;

    extends DynamicVCC.Components.Units.HX.FrostDefrostHX.BaseClasses.PartialFrostGrowth;

    parameter SI.Density rho_ice=917;

    SI.Density rho_wf "Water vapor density at frost surface";
    SI.Density rho_wp "Water vapor density at tube surface";
    SI.DiffusionCoefficient D(min=0) "Effective diffusion coefficient";
    SI.DiffusionCoefficient D_s(min=0) "Diffusion coefficient water-vapor";
    Real epsilon(min=0.0,max=1.0) "Porosity";
    Real tau "Tortuosity";

    output Real m_rho(min=0);
    output Real m_x(min=0);

  protected
    Real gamma(min=0);
    Real alpha_f(min=0);
    Real m_f "Total mass flux";

    //SI.Thickness x_f_nominal(start=x_f_init);
    //SI.Density rho_f_nominal(start=rho_f_init);

  equation

    // Properties
    rho_wf=waterVaporDensity_sat(T_fs);
    rho_wp=waterVaporDensity_sat(T_ps);
    k_f=thermalConductivity(rho_f);
    D_s=2.302*0.98e5/1e5*(T_fs/256)^1.81/1e5;
    epsilon=1-rho_f/rho_ice;
    tau=epsilon/(1-(1-epsilon)^0.5);
    D=D_s*epsilon/tau;

    // heat and mass transfer
    gamma=acosh(max(1+1e-6,rho_wf/rho_wp))/x_f;
    alpha_f=gamma^2*D;
    T_fs=T_ps+q_dot_tot/k_f*x_f-alpha_f*delta_h_ig*rho_wp*(cosh(gamma*x_f)-1)/(k_f*gamma^2);

    m_f=m_flow_dehumid/As;
    m_rho=max(0,min(m_f,D*gamma*rho_wf*tanh(gamma*x_f)));
    m_x=m_f-m_rho;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end FrostGrowth_Lee;

  model FrostMelt_Fuzzy "1-D frost melting model transitioned using Fuzzy logic"

    import Modelica.Constants.eps;
    import Modelica.Fluid.Utilities.regStep;
    import DynamicVCC.Components.Units.HX.FrostDefrostHX.Utilities.thermalConductivity_Lee;

    replaceable package MediumWater=Modelica.Media.Water.WaterIF97_pT;

    replaceable package MediumAir=Modelica.Media.Air.MoistAir;

  /************ Connectors *******************/
    DynamicVCC.Interfaces.HeatPort_in heatPort_a;
    DynamicVCC.Interfaces.HeatPort_in heatPort_b;

    input SI.Area As "Surface area";
    input SI.Thickness tau_t "Tube wall thickness";
    input SI.Temperature T_amb;
    input Real x_steam;

    parameter SI.Thickness tau_small=1e-7;

  /************* Properties ********************/

    parameter SI.Temperature T_melt=273.15;
    parameter SI.Density rho_ice=917, rho_water=997;
    parameter SI.SpecificHeatCapacity cp_water=4190;
    parameter SI.ThermalConductivity k_water=0.5826;
    parameter SI.SpecificHeatCapacity cp_ice=2126;
    parameter SI.SpecificEnergy L_fus=333606 "Water heat of fusion";
    parameter SI.SpecificEnergy L_fg=2.497e6 "Water heat of vaporization";
    parameter SI.Thickness tau_water_max=0.02e-3 "Maximum water thickness held";
    parameter SI.ThermalConductivity k_t=386 "Tube wall conductivity";
    parameter SI.Pressure p_atm=MediumAir.p_default;

    // Initialization
    parameter SI.Temperature T_f_init=MediumWater.reference_T;
    parameter SI.Thickness tau_f_init=1e-3;
    parameter SI.Density rho_f_init=300;

  /************* Variables ********************/
    MediumAir.ThermodynamicState airState;
    MediumAir.Temperature T_air(start=MediumAir.T_default);
    MediumAir.Density rho_air;
    Real phi_air "Relative humidity";
    MediumAir.SpecificHeatCapacity cp_air;
    MediumAir.ThermalConductivity k_air;
    MediumAir.MassFraction X_air[MediumAir.nX]={x_steam,1-x_steam};
    SI.Thickness tau_air(start=tau_small,fixed=true);
    SI.Pressure pv "Partial pressure of water vapor";
    SI.Density rhov "Density of water vapor";
    SI.Density rhov_s "Saturation density of water vapor at water film surface temperature";

    SI.Temperature T_f(start=T_f_init,fixed=true);
    SI.Temperature T_fs(start=T_f_init) "Frost surface temperature";
    input SI.Thickness tau_f "Frost thickness";
    SI.ThermalConductivity k_f;
    SI.SpecificHeatCapacity cp_f;
    input SI.Density rho_f;
    Real epsilon=1-rho_f/rho_ice "porosity";

    SI.Temperature T_water(min=T_melt,start=T_melt,fixed=true);
    SI.Temperature T_water_s(start=T_melt) "Water film surface temperature in vaporizing stage";
    SI.Thickness tau_water(start=tau_small,fixed=true);

    SI.Temperature T_t "Tube or fin wall temperature";
    SI.Temperature T_tf "Tube-frost surface temperature";
    SI.Temperature T_tw "Tube-water surface temperature";
    SI.EnergyFlowRate Q_t "Heat transfer rate through tube wall";
    SI.Temperature T_as "Temperature of surface with ambient";

  protected
    Real xdot_preheat[6];
    Real xdot_meltstart[6];
    Real xdot_melt[6];
    Real xdot_vapor[6];
    Real xdot_dryheat[6];
    Real xdot_fuzzy[6]; //T_f,tau_f,T_water,tau_water,T_a,tau_a

  /************************ Fuzzy Model ***************************/
    Real mu_Tw_N; //Membership function
    Real mu_Tw_P;
    Real mu_tauWater_N;
    Real mu_tauWater_P;
    Real mu_tauWater_LP;
    Real mu_tauf_N;
    Real mu_tauf_P;

    Real mu_rule[5]; //Fuzzy rules
    Real weights[5];

    Real Q_water; //Water film exist
    Real Q_preheat;

  equation
  /************ Properties calculation *******************/

    rho_air=MediumAir.density(airState);
    airState=MediumAir.setState_pTX(p_atm,T_amb,X_air);
    phi_air=MediumAir.relativeHumidity(airState);
    cp_air=MediumAir.specificHeatCapacityCp(airState);
    k_air=MediumAir.thermalConductivity(airState);
    pv=phi_air*MediumAir.saturationPressure(airState.T);
    rhov=pv/(461.5*T_amb) "Ideal gas law for water vapor";
    rhov_s=2.55e-2;//Air.saturationPressure(max(T_water,T_melt))/(461.5*T_water);

    //waterState=Water.setState_pT(Water.p_default,max(T_melt,T_water));
    //cp_water=Water.specificHeatCapacityCp(waterState);
    //k_water=Water.thermalConductivity(waterState);

    //rho_f=rho_f_init "Assume constant frost density throughout melting";
    k_f=thermalConductivity_Lee(rho_f);
    cp_f=(epsilon*rho_air*cp_air+(1-epsilon)*rho_ice*cp_ice)/rho_f;

    //preheating
    k_t*(T_t-T_tf)/(tau_t/2)=k_f*(T_tf-T_f)/(tau_f/2);
    rho_f*cp_f*tau_f*xdot_preheat[1] = k_f*(T_tf+T_fs-2*T_f)/(tau_f/2);
    xdot_preheat[2:end]=zeros(5);
    //k_f*(T_f-T_fs)/(tau_f/2)=alpha*(T_fs-T_amb);
    k_f*(T_f-T_fs)/(tau_f/2)=-heatPort_b.Q_flow/As;

    //melting start
    k_t*(T_t-T_tw)/(tau_t/2)=k_water*(T_tw-T_water)/(tau_water/2);
    rho_f*cp_f*tau_f*xdot_meltstart[1]=k_f*(T_melt+T_fs-2*T_f)/(tau_f/2);
    rho_water*cp_water*tau_water*xdot_meltstart[3] = k_water*(T_melt+T_tw-2*T_water)/(tau_water/2);
    -rho_f*L_fus*xdot_meltstart[2]=k_water*(T_water-T_melt)/(tau_water/2)-k_f*(T_melt-T_f)/(tau_f/2);
    xdot_meltstart[4]=-xdot_meltstart[2];
    xdot_meltstart[5:end]=zeros(2);

    //melting
    rho_f*cp_f*tau_f*xdot_melt[1]=k_f*(T_melt+T_fs-2*T_f)/(tau_f/2);
    rho_water*cp_water*tau_water*xdot_melt[3] = k_water*(T_tw-T_water)/(tau_water/2) - (T_water-T_air)/(tau_air/2/k_air+tau_water/2/k_water);
    rho_air*cp_air*tau_air*xdot_melt[5] = (T_water-T_air)/(tau_air/2/k_air+tau_water/2/k_water) - k_air*(T_air-T_melt)/(tau_air/2);
    -rho_f*L_fus*xdot_melt[2]=k_air*(T_air-T_melt)/(tau_air/2)-k_f*(T_melt-T_f)/(tau_f/2);
    xdot_melt[6] = -xdot_melt[2];
    xdot_melt[4]=0;

    //vaporizing

    xdot_vapor[1:2]=zeros(2);
    rho_water*cp_water*tau_water*xdot_vapor[3] = k_water*(T_tw+T_water_s-2*T_water)/(tau_water/2);
    rho_water*xdot_vapor[4] = 0.0085*(rhov-rhov_s); //mass transfer
    //k_water*(T_water-T_water_s)/(tau_water/2)=alpha*(T_water_s-T_amb)+0.0085*(rhov-rhov_s)*L_fg;
    k_water*(T_water-T_water_s)/(tau_water/2)=-heatPort_b.Q_flow/As+0.0085*(rhov-rhov_s)*L_fg;
    xdot_vapor[5:end]=zeros(2);

    //dry heating
    xdot_dryheat=zeros(6);

  /************************** Fuzzy modeling *************************/
    //mu_Tw_N=regStep(T_w-T_melt,0,1,T_melt*1e-3);
    //mu_Tw_P=regStep(T_w-T_melt,1,0,T_melt*1e-3);
    //mu_Tw_N=smooth(0,noEvent(if T_w<T_melt then 1 elseif T_w>T_melt+T_melt*1e-5 then 0 else 1-(T_w-T_melt)/(T_melt*1e-5)));
    mu_Tw_P=smooth(0,noEvent(if T_t<T_melt then 0 elseif T_t>T_melt+T_melt*1e-4 then 1 else (T_t-T_melt)/(T_melt*1e-4)));
    //mu_Tw_P=smooth(0,noEvent(if T_w<T_melt then 0 elseif T_w>T_melt+T_melt*1e-5 then 1 else ((T_w-T_melt)/(T_melt*1e-5))^3));
    mu_Tw_N=smooth(0,noEvent(if T_t<T_melt then 1 elseif T_t>T_melt+T_melt*1e-4 then 0 else (-(T_t-T_melt-T_melt*1e-4)/(T_melt*1e-4))^5));

    //mu_tauWater_N=regStep(tau_water-tau_water_max*1e-4,0,1,tau_water_max*1e-4);
    mu_tauWater_N=smooth(0,noEvent(if tau_water<tau_water_max*1e-6 then 1 elseif tau_water>2*tau_water_max*1e-6 then 0 else 1-(tau_water-tau_water_max*1e-6)/(tau_water_max*1e-6)));

    //mu_tauWater_LP=regStep(tau_water_max-tau_water,0,1,tau_water_max*1e-4);
    mu_tauWater_LP=smooth(0,noEvent(if tau_water>tau_water_max then 1 elseif tau_water<tau_water_max-tau_water_max*1e-6 then 0 else (tau_water-tau_water_max*(1-1e-6))/(tau_water_max*1e-6)));

  /*
  if tau_water>tau_water_max/2 then
    //mu_tauWater_P=regStep(tau_water_max-tau_water,1,0,tau_water_max*1e-4);
    mu_tauWater_P=smooth(0,noEvent(if tau_water<tau_water_max-tau_water_max*1e-5 then 1 elseif tau_water>tau_water_max then 0 else 1-(tau_water-tau_water_max*(1-1e-5))/(tau_water_max*1e-5)));
  else
    //mu_tauWater_P=regStep(tau_water-tau_water_max*1e-4,1,0,tau_water_max*1e-4);
    mu_tauWater_P=smooth(0,noEvent(if tau_water<eps then 0 elseif tau_water>tau_water_max*1e-4 then 1 else tau_water/(tau_water_max*1e-4)));
    end if;
    */

      mu_tauWater_P=smooth(0,noEvent(if tau_water<eps then 0 elseif tau_water<tau_water_max*1e-6 then tau_water/(tau_water_max*1e-6) elseif tau_water<tau_water_max-tau_water_max*1e-6 then 1 elseif tau_water<tau_water_max then
      1-(tau_water-tau_water_max*(1-1e-6))/(tau_water_max*1e-6) else 0));

    //mu_tauf_N=regStep(tau_f-tau_water_max*1e-6,0,1,tau_water_max*1e-4);
    //mu_tauf_P=regStep(tau_f-tau_water_max*1e-6,1,0,tau_water_max*1e-4);
    mu_tauf_N=smooth(0,noEvent(if tau_f<1e-8 then 1 elseif tau_f>2e-8 then 0 else 1-(tau_f-1e-8)/1e-8));
    mu_tauf_P=smooth(0,noEvent(if tau_f<1e-8 then 0 elseif tau_f>2e-8 then 1 else (tau_f-1e-8)/1e-8));

    mu_rule[1]=mu_Tw_N;
    mu_rule[2]=mu_Tw_P*(mu_tauWater_N+mu_tauWater_P)*mu_tauf_P;
    mu_rule[3]=mu_Tw_P*mu_tauWater_LP*mu_tauf_P;
    mu_rule[4]=mu_Tw_P*(mu_tauWater_P+mu_tauWater_LP)*mu_tauf_N;
    mu_rule[5]=mu_Tw_P*mu_tauWater_N*mu_tauf_N;

    weights=mu_rule/sum(mu_rule);

    xdot_fuzzy=xdot_preheat.*weights[1]+xdot_meltstart.*weights[2]+xdot_melt.*weights[3]+xdot_vapor.*weights[4]+xdot_dryheat.*weights[5];
    der(T_f)=xdot_fuzzy[1];
    //der(tau_f)=xdot_fuzzy[2];
    der(T_water)=xdot_fuzzy[3];
    der(tau_water)=xdot_fuzzy[4];
    der(T_air)=xdot_fuzzy[5];
    der(tau_air)=xdot_fuzzy[6];

    Q_water=k_water*As*(T_t-T_water)/(tau_water/2);
    Q_preheat=k_f*As*(T_tf-T_f)/(tau_f/2);
    //Q_t=-heatPort_b.Q_flow*(weights[1]+weights[5])+Q_water*sum(weights[2:4]);
    T_as=T_fs*sum(weights[1:3])+T_water_s*weights[4]+T_t*weights[5];
    Q_t=-heatPort_b.Q_flow;
    //k_water*(T_water-T_water_s_vaporize)/(tau_water/2) = HTC*(T_water_s_vaporize-T_amb) - rho_water*L_fg*der(tau_water);
    //T_water_s=T_water*sum(weights[1:3])+T_water_s_vaporize*weights[4]+T_water*weights[5];

  /************************** Boundary conditions ******************************/
    T_t=heatPort_a.T;
    Q_t=heatPort_a.Q_flow;
    heatPort_b.T=T_as;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=2000, __Dymola_Algorithm="Dassl"));
  end FrostMelt_Fuzzy;

  model FrostMelt_Fuzzy_new "1-D frost melting model transitioned using Fuzzy logic 01202022"

    import Modelica.Constants.eps;
    import Modelica.Fluid.Utilities.regStep;
    import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
    import thermalConductivity=DynamicVCC.Components.Units.HX.FrostDefrostHX.Utilities.thermalConductivity_Lee;
    import DynamicVCC.Components.Units.HX.FrostDefrostHX.Utilities.waterVaporDensity_sat;

    input SI.Temperature T_t;

    input SI.Temperature T_amb;

    output SI.EnergyFlowRate Q_flow "Heat transfer rate through tube wall";

    output Real xdot_fuzzy[7]; //T_f,tau_f,T_water,tau_water,T_a,tau_a,rho_f

    replaceable package MediumWater=Modelica.Media.Water.WaterIF97_pT;

    replaceable package MediumAir=Modelica.Media.Air.MoistAir;

  /************ Connectors *******************/
    //TransientVCC.Interfaces.HeatPort_in heatPort_a;
    //TransientVCC.Interfaces.HeatPort_in heatPort_b;

    input SI.Area As "Surface area";
    input MediumAir.ThermodynamicState airState;
    input SI.Thickness tau_f(start=tau_f_init) "Frost thickness";
    input SI.Density rho_f(start=rho_f_init) "Frost density";

    parameter SI.Thickness tau_small=1e-5;

  /************* Properties ********************/

    parameter SI.Temperature T_melt=273.15;
    parameter SI.Density rho_ice=917, rho_water=997;
    parameter SI.SpecificHeatCapacity cp_water=4190;
    parameter SI.ThermalConductivity k_water=0.5826;
    parameter SI.SpecificHeatCapacity cp_ice=2126;
    parameter SI.SpecificEnergy L_fus=333606 "Water heat of fusion";
    parameter SI.SpecificEnergy L_fg=2.497e6 "Water heat of vaporization";
    parameter SI.Thickness tau_water_max=0.015e-3 "Maximum water thickness held";
    parameter Real c_v=0.05 "Water vaporization coefficient";
    parameter SI.ThermalConductivity k_t=386 "Tube wall conductivity";
    parameter SI.Pressure p_atm=MediumAir.p_default;
    parameter SI.CoefficientOfHeatTransfer h_free=0.5 "Free convection";

    // Initialization
    parameter SI.Temperature T_f_init=T_melt;
    parameter SI.Thickness tau_f_init=1e-3;
    parameter SI.Density rho_f_init=300;

  /************* Variables ********************/
    MediumAir.Temperature T_air(start=T_f_init);
    MediumAir.Density rho_air;
    Real phi_air "Relative humidity";
    MediumAir.SpecificHeatCapacity cp_air;
    MediumAir.ThermalConductivity k_air;
    SI.Thickness tau_air(start=tau_small,fixed=true);
    SI.Pressure pv "Partial pressure of water vapor";
    SI.Density rhov "Density of water vapor";
    SI.Density rhov_s "Saturation density of water vapor at water film surface temperature";

    SI.Temperature T_f(start=T_f_init,fixed=true);
    SI.Temperature T_fs(start=T_f_init) "Frost surface temperature";
    SI.ThermalConductivity k_f;
    SI.SpecificHeatCapacity cp_f;

    Real epsilon=1-rho_f/rho_ice "porosity";

    SI.Temperature T_water(start=T_f_init,fixed=true);
    SI.Temperature T_water_s(start=T_melt) "Water film surface temperature in vaporizing stage";
    SI.Thickness tau_water(start=tau_small,fixed=true);
    SI.Temperature T_as "Temperature of surface with ambient";

  protected
    Real xdot_preheat[7];
    Real xdot_meltstart[7];
    Real xdot_melt[7];
    Real xdot_vapor[7];
    Real xdot_dryheat[7];

  /************************ Fuzzy Model ***************************/
    Real mu_Tw_N; //Membership function
    Real mu_Tw_P;
    Real mu_tauWater_N;
    Real mu_tauWater_P;
    Real mu_tauWater_LP;
    Real mu_tauf_N;
    Real mu_tauf_P;

    Real mu_rule[5]; //Fuzzy rules
    Real weights[5];

    Real Q_water; //Water film exist
    Real Q_preheat;

  equation
  /************ Properties calculation *******************/

    rho_air=MediumAir.density(airState);
    phi_air=MediumAir.relativeHumidity(airState);
    cp_air=MediumAir.specificHeatCapacityCp(airState);
    k_air=MediumAir.thermalConductivity(airState);
    pv=phi_air*MediumAir.saturationPressure(MediumAir.temperature(airState));
    rhov=pv/(461.5*MediumAir.temperature(airState)) "Ideal gas law for water vapor";
    rhov_s=waterVaporDensity_sat(T_as);

    //waterState=Water.setState_pT(Water.p_default,max(T_melt,T_water));
    //cp_water=Water.specificHeatCapacityCp(waterState);
    //k_water=Water.thermalConductivity(waterState);

    //rho_f=rho_f_init "Assume constant frost density throughout melting";
    k_f=thermalConductivity(rho_f);
    cp_f=(epsilon*rho_air*cp_air+(1-epsilon)*rho_ice*cp_ice)/rho_f;

    //preheating
    rho_f*cp_f*tau_f*xdot_preheat[1] = k_f*(T_t+T_fs-2*T_f)/(tau_f/2);
    xdot_preheat[2]=0;
    xdot_preheat[3]=(T_t-T_water)/1e-5;
    xdot_preheat[4]=0;
    xdot_preheat[5]=(T_t-T_air)/1e-4;
    xdot_preheat[6:7]=zeros(2);
    k_f*(T_f-T_fs)/(tau_f/2)=h_free*(T_fs-T_amb);
    //1e3*(T_f-T_fs)/(tau_f/2)=-heatPort_b.Q_flow/As;

    //melting start
    rho_f*cp_f*tau_f*xdot_meltstart[1]=k_f*(T_melt+T_fs-2*T_f)/(tau_f/2);
    rho_water*cp_water*tau_water*xdot_meltstart[3] = k_water*(T_melt+T_t-2*T_water)/(tau_water/2);
    -rho_f*L_fus*xdot_meltstart[2]=k_water*(T_water-T_melt)/(tau_water/2)-k_f*max(0,(T_melt-T_f))/(tau_f/2);
    xdot_meltstart[4]*rho_water=-xdot_meltstart[2]*rho_f; // Mass balance at interface
    xdot_meltstart[5]=(T_t-T_air)/1e-4;
    xdot_meltstart[6]=0;
    xdot_meltstart[7]=0;

    //melting
    rho_f*cp_f*tau_f*xdot_melt[1]=k_f*(T_melt+T_fs-2*T_f)/(tau_f/2);
    rho_water*cp_water*tau_water*xdot_melt[3] = k_water*(T_t-T_water)/(tau_water/2) - (T_water-T_air)/(tau_air/2/k_air+tau_water/2/k_water);
    rho_air*cp_air*tau_air*xdot_melt[5] = (T_water-T_air)/(tau_air/2/k_air+tau_water/2/k_water) - k_air*(T_air-T_melt)/(tau_air/2);
    -rho_f*L_fus*xdot_melt[2]=k_air*(T_air-T_melt)/(tau_air/2)-k_f*(T_melt-T_f)/(tau_f/2);
    xdot_melt[6] = -xdot_melt[2];
    xdot_melt[4]=0;
    xdot_melt[7]=0;

    //vaporizing

    xdot_vapor[1:2]=zeros(2);
    rho_water*cp_water*tau_water*xdot_vapor[3] = k_water*(T_t+T_water_s-2*T_water)/(tau_water/2);
    rho_water*xdot_vapor[4] = c_v*(rhov-rhov_s); //mass transfer
    k_water*(T_water-T_water_s)/(tau_water/2)=h_free*(T_water_s-T_amb)+c_v*(rhov_s-rhov)*L_fg;
    xdot_vapor[5:6]=zeros(2);
    xdot_vapor[7]=(25-rho_f)/1e-2;
    xdot_dryheat[1:6]=zeros(6);

    //dry heating
    xdot_dryheat[7]=(25-rho_f)/1e-2;

  /************************** Fuzzy modeling *************************/
    //mu_Tw_N=regStep(T_w-T_melt,0,1,T_melt*1e-3);
    //mu_Tw_P=regStep(T_w-T_melt,1,0,T_melt*1e-3);
    //mu_Tw_N=smooth(0,noEvent(if T_w<T_melt then 1 elseif T_w>T_melt+T_melt*1e-5 then 0 else 1-(T_w-T_melt)/(T_melt*1e-5)));
    //mu_Tw_P=smooth(0,noEvent(if T_t<T_melt then 0 elseif T_t>T_melt+1e-2 then 1 else (T_t-T_melt)/(T_melt*1e-4)));
    mu_Tw_P=spliceFunction(1.0,0.0,T_t-T_melt-1e-2,1e-2);
    //mu_Tw_N=smooth(0,noEvent(if T_t<T_melt then 1 elseif T_t>T_melt+T_melt*1e-4 then 0 else (-(T_t-T_melt-T_melt*1e-4)/(T_melt*1e-4))^5));
    mu_Tw_N=spliceFunction(0.0,1.0,T_t-T_melt-1e-2,1e-2);

    //mu_tauWater_N=regStep(tau_water-tau_water_max*1e-4,0,1,tau_water_max*1e-4);
    //mu_tauWater_N=smooth(0,noEvent(if tau_water<tau_water_max*1e-2 then 1 elseif tau_water>2*tau_water_max*1e-2 then 0 else 1-(tau_water-tau_water_max*1e-2)/(tau_water_max*1e-2)));
    mu_tauWater_N=smooth(0,noEvent(if tau_water<1*1e-5 then 1 elseif tau_water>1.1e-5 then 0 else 1-(tau_water-1e-5)/(1e-6)));

    //mu_tauWater_LP=regStep(tau_water_max-tau_water,0,1,tau_water_max*1e-4);
    //mu_tauWater_LP=smooth(0,noEvent(if tau_water>tau_water_max then 1 elseif tau_water<tau_water_max-tau_water_max*1e-6 then 0 else (tau_water-tau_water_max*(1-1e-6))/(tau_water_max*1e-6)));
    //mu_tauWater_LP=smooth(0,noEvent(if tau_water<tau_water_max*(1-1e-2) then 0 elseif tau_water<tau_water_max then (tau_water-tau_water_max*(1-1e-2))/(tau_water_max*1e-2) else 1));
    mu_tauWater_LP=smooth(0,noEvent(if tau_water<tau_water_max*(1-5e-2) then 0 elseif tau_water<tau_water_max then (tau_water-tau_water_max*(1-5e-2))/(tau_water_max*5e-2) else 1));

  /*
  if tau_water>tau_water_max/2 then
    //mu_tauWater_P=regStep(tau_water_max-tau_water,1,0,tau_water_max*1e-4);
    mu_tauWater_P=smooth(0,noEvent(if tau_water<tau_water_max-tau_water_max*1e-5 then 1 elseif tau_water>tau_water_max then 0 else 1-(tau_water-tau_water_max*(1-1e-5))/(tau_water_max*1e-5)));
  else
    //mu_tauWater_P=regStep(tau_water-tau_water_max*1e-4,1,0,tau_water_max*1e-4);
    mu_tauWater_P=smooth(0,noEvent(if tau_water<eps then 0 elseif tau_water>tau_water_max*1e-4 then 1 else tau_water/(tau_water_max*1e-4)));
    end if;
    */

    //mu_tauWater_P=smooth(0,noEvent(if tau_water<tau_water_max*1e-2 then 0 elseif tau_water<2*tau_water_max*1e-2 then (tau_water-tau_water_max*1e-2)/(tau_water_max*1e-2) elseif tau_water<tau_water_max*(1-1e-2) then 1 elseif tau_water<tau_water_max then
      //1-(tau_water-tau_water_max*(1-1e-2))/(tau_water_max*1e-2) else 0));
    mu_tauWater_P=smooth(0,noEvent(if tau_water<1*1e-5 then 0 elseif tau_water<1.1e-5 then (tau_water-1*1e-5)/(1e-6) elseif tau_water<tau_water_max*(1-5e-2) then 1 elseif tau_water<tau_water_max then
      1-(tau_water-tau_water_max*(1-5e-2))/(tau_water_max*5e-2) else 0));

    //mu_tauf_N=regStep(tau_f-tau_water_max*1e-6,0,1,tau_water_max*1e-4);
    //mu_tauf_P=regStep(tau_f-tau_water_max*1e-6,1,0,tau_water_max*1e-4);
    //mu_tauf_N=smooth(0,noEvent(if tau_f<tau_small then 1 elseif tau_f>2*tau_small then 0 else 1-(tau_f-tau_small)/tau_small));
    //mu_tauf_P=smooth(0,noEvent(if tau_f<tau_small then 0 elseif tau_f>2*tau_small then 1 else (tau_f-tau_small)/tau_small));
    mu_tauf_N=smooth(0,noEvent(if tau_f<1e-5 then 1 elseif tau_f>1.1*1e-5 then 0 else 1-(tau_f-1e-5)/1e-6));
    mu_tauf_P=smooth(0,noEvent(if tau_f<1e-5 then 0 elseif tau_f>1.1*1e-5 then 1 else (tau_f-1e-5)/1e-6));

    mu_rule[1]=mu_Tw_N;
    mu_rule[2]=mu_Tw_P*(mu_tauWater_N+mu_tauWater_P)*mu_tauf_P;
    mu_rule[3]=mu_Tw_P*mu_tauWater_LP*mu_tauf_P;
    mu_rule[4]=mu_Tw_P*(mu_tauWater_P+mu_tauWater_LP)*mu_tauf_N;
    mu_rule[5]=mu_Tw_P*mu_tauWater_N*mu_tauf_N;

    weights=mu_rule/sum(mu_rule);

    xdot_fuzzy=xdot_preheat.*weights[1]+xdot_meltstart.*weights[2]+xdot_melt.*weights[3]+xdot_vapor.*weights[4]+xdot_dryheat.*weights[5];
    der(T_f)=xdot_fuzzy[1];
    der(T_water)=xdot_fuzzy[3];
    der(tau_water)=xdot_fuzzy[4];
    der(T_air)=xdot_fuzzy[5];
    der(tau_air)=xdot_fuzzy[6];

    Q_water=k_water*As*(T_t-T_water)/(tau_water/2);
    Q_preheat=k_f*As*(T_t-T_fs)/tau_f;
    Q_flow=As*h_free*(T_t-T_amb)*(weights[1]+weights[5])+Q_water*sum(weights[2:4]);
    T_as=T_fs*sum(weights[1:3])+T_water_s*weights[4]+T_t*weights[5];
    //Q_t=-heatPort_b.Q_flow;
    //k_water*(T_water-T_water_s_vaporize)/(tau_water/2) = HTC*(T_water_s_vaporize-T_amb) - rho_water*L_fg*der(tau_water);
    //T_water_s=T_water*sum(weights[1:3])+T_water_s_vaporize*weights[4]+T_water*weights[5];

  /************************** Boundary conditions ******************************/
    //T_t=heatPort_a.T;
    //Q_t=heatPort_a.Q_flow;
    //heatPort_b.T=T_as;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=2000, __Dymola_Algorithm="Dassl"));
  end FrostMelt_Fuzzy_new;

  package BaseClasses "Base classes of 1-D frost growth models"
    extends Modelica.Icons.BasesPackage;
    partial model PartialFrostGrowth

      import Modelica.Math.Vectors.interpolate;
      import Modelica.Fluid.Utilities.regStep;

      package Medium=Modelica.Media.Air.MoistAir;

      input SI.Area As "Surface Area";

      input Medium.MassFlowRate m_flow_dehumid "Mass transfer rate to frost";

      input SI.Thickness x_f(start=x_f_init) "Frost thickness";

      input SI.Density rho_f(start=rho_f_init);

      input SI.Temperature T_t;

      output SI.HeatFlowRate Q_flow;

      parameter SI.Thickness x_f_init=1e-5;
      parameter SI.Density rho_f_init=30;
      parameter SI.Time T_sample=10;
      parameter SI.SpecificEnthalpy delta_h_ig=2836.6e+3 "Latent heat of sublimation, approximated to a constant";

      // Connectors
      DynamicVCC.Interfaces.HeatPort_out heatPort "Heat transfer to air";

      Medium.Temperature T_ps "Cold plate surface temperature";
      Medium.Temperature T_fs "Frost surface temperature";

    protected
      SI.ThermalConductivity k_f "Frost layer thermal conductivity is a function of the frost density";
      SI.HeatFlux q_dot_tot "Total heat flux on frost surface";
      //SI.Temperature T_ps_in(start=274);
    equation
      //der(T_ps)=(min(275,T_ps_in)-T_ps)/1e-5;
      T_ps=T_t;
      // Boundary conditions
      //heatPort_a.T=T_ps;
      heatPort.T=T_fs;
      q_dot_tot=heatPort.Q_flow/As;
      //heatPort_a.Q_flow=-q_dot_tot*As "Neglect energy storage of frost layer";
      //heatPort_a.Q_flow=-k_f*As*(T_fs-T_ps)/x_f "Consider thermal resistence of the frost layer";
      Q_flow=-k_f*As*(T_fs-T_ps)/x_f;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialFrostGrowth;
  end BaseClasses;

  package Utilities
    extends Modelica.Icons.UtilitiesPackage;
    function thermalConductivity_Lee "Frost thermal conductivity correlation (Lee, 1997)"
      extends Modelica.Icons.Function;
      input Real rho;
      output SI.ThermalConductivity lambda;
    algorithm
      lambda:=0.132 + 3.13e-4*rho + 1.6e-7*rho^2;
    end thermalConductivity_Lee;

    function waterVaporDensity_sat "Saturated water vapor density as a function of temperature between 250 and 300K"
      import Modelica.Media.Air.MoistAir.Utilities.spliceFunction;
      input SI.Temperature T;
      output SI.Density rho;
    protected
      SI.Density rho_sub "sublimation line";
      SI.Density rho_vapor "vapor-liquid line";
    algorithm
      rho_sub:=6.539e-06*T^2 - 0.003247*T + 0.4038;
      rho_vapor:=2.008e-05*T^2 - 0.01076*T + 1.445;
      rho:=spliceFunction(
              rho_vapor,
              rho_sub,
              T - 273.15,
              0.1);
    end waterVaporDensity_sat;
  end Utilities;
end FrostDefrostHX;
