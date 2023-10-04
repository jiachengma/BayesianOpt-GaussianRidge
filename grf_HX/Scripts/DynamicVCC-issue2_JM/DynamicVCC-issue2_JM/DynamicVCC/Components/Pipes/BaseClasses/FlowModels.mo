within DynamicVCC.Components.Pipes.BaseClasses;
package FlowModels "Flow models for evaluating momentum balances"
  extends Modelica.Icons.VariantsPackage;
  partial model PartialStaggeredGridFlowModel "Partial staggered grid flow models (momentum grid)"

    import Modelica.Fluid.Utilities.regStep;

    // extendisng flow cells
    extends DynamicVCC.Components.Pipes.BaseClasses.PartialFlowCell(final m=n - 1);

    parameter Integer n=2 "Number of control volume";

    // Inputs to flow cells
    input Medium.ThermodynamicState states[n];

    input SI.Velocity vs[n] "Flow velocity at cell boundaries";

    input SI.Area crossAreas[n] "Cross flow area";

    input SI.Length dimensions[n] "Characteristic dimension of fluid flow, typically pipe diameter";

    // Numerical settings
    parameter Boolean use_I_flows= momentumDynamics<>DynamicVCC.Components.Types.Dynamics.SteadyState
    "= true to evaluate momentum flow differences in momentum balances";
    parameter Boolean upStreamProperties=true
    "= true to use upstream thermodynamic properties, otherwise use average properties";
    parameter Boolean from_dp=momentumDynamics==DynamicVCC.Components.Types.Dynamics.SteadyState
    "= true, compute m_flow from dp, otherwise compute dp from m_flow";


    // Initial conditions
    parameter Medium.AbsolutePressure p_a_start=Medium.p_default
    "Start value for pressure at design inflow";
    parameter Medium.AbsolutePressure p_b_start=Medium.p_default
    "Start value for pressure at design outflow";

    // Variables
    Medium.Density rho[n]=Medium.density(states);
    Medium.Density rho_act[n-1] "crossAreastual density of flow cell";
    Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states);
    Medium.DynamicViscosity mu_act[n-1] "crossAreastual viscosity of flow cell";

    SI.Pressure dps_f[n-1](each start=(p_a_start-p_b_start)/(n-1)) "Frictional pressure drop";

    // Parameters
    parameter SI.AbsolutePressure dp_nominal=1e3*dp_small
    "Nominal pressure drop";
    parameter SI.MassFlowRate m_flow_nominal=system.m_flow_nominal
    "Nominal mass flow rate";
    parameter SI.MassFlowRate m_flow_small=system.m_flow_small;

    parameter Boolean use_dp_nominal=false
    "Initial dp_nominal from flow model";

    parameter SI.AbsolutePressure dp_small=system.dp_small;

    SI.Diameter diameters[n-1]=0.5*(dimensions[1:n-1]+dimensions[2:n])
    "Mean diameters of flow cells";


  protected
    parameter Medium.Density rho_nominal=Medium.density(
      Medium.setState_phX(Medium.p_default,Medium.h_default));
    parameter Medium.DynamicViscosity mu_nominal=Medium.dynamicViscosity(
      Medium.setState_phX(Medium.p_default,Medium.h_default));
  equation

    // Properties
    if not enableReverseFlow then
      rho_act=rho[1:n-1];
      mu_act=mu[1:n-1];
    elseif not upStreamProperties then
      rho_act=0.5*(rho[1:n-1]+rho[2:n]);
      mu_act=0.5*(mu[1:n-1]+mu[2:n]);
    else
      for i in 1:n-1 loop
        rho_act[i]=noEvent(if m_flows[i]>0 then rho[i] else rho[i+1]);
        mu_act[i]=noEvent(if m_flows[i]>0 then mu[i] else mu[i+1]);
      end for;
    end if;

    // Source terms of momentum balances
    if use_I_flows then
      I_flows={rho[i]*vs[i]^2*crossAreas[i] - rho[i+1]*vs[i+1]^2*crossAreas[i+1] for i in 1:n-1};
    else
      I_flows=zeros(n-1);
    end if;

    F_p={0.5*(crossAreas[i]+crossAreas[i+1])*(Medium.pressure(states[i+1])-Medium.pressure(states[i])) for i in 1:n-1};

    dps_f={F_f[i]*2/(crossAreas[i]+crossAreas[i+1]) for i in 1:n-1};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                            Line(
            points={{-80,-60},{-80,60},{80,-60},{80,62}},
            color={0,0,255},
            thickness=1)}),                                        Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialStaggeredGridFlowModel;

  model DetailedFlow "Staggered grid for solving momentum balances with detailed friction model"

    extends DynamicVCC.Components.Pipes.BaseClasses.FlowModels.PartialStaggeredGridFlowModel;

    replaceable package WallFriction=DynamicVCC.Components.Pipes.BaseClasses.WallFriction.Correlation_SinglePhase
      constrainedby DynamicVCC.Components.Pipes.BaseClasses.WallFriction.PartialWallFriction;


  protected
    SI.AbsolutePressure dp_f_nominal=
    sum(WallFriction.pressureLoss_m_flow(
      m_flow_nominal,
      rho_nominal,
      mu_nominal,
      lengths,
      diameters,
      (crossAreas[1:n-1]+crossAreas[2:n])/2,
      m_flow_small))
      "Nominal pressure drop";

  equation
    if from_dp then
      m_flows=homotopy(
        actual=WallFriction.massFlowRate_dp(
        dps_f,
        rho_act,
        mu_act,
        lengths,
        diameters,
        (crossAreas[1:n-1]+crossAreas[2:n])/2,
        dp_small),
        simplified=dps_f/dp_nominal*m_flow_nominal);
    else
      dps_f=homotopy(
        actual=WallFriction.pressureLoss_m_flow(
        m_flows,
        rho_act,
        mu_act,
        lengths,
        diameters,
        (crossAreas[1:n-1]+crossAreas[2:n])/2,
        m_flow_small),
        simplified=m_flows/m_flow_nominal*dp_nominal);
    end if;
  initial equation
    //dp_nominal=dp_f_nominal;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end DetailedFlow;

  model ConstantFrictionFlow "Constant friction factor"

    import Modelica.Fluid.Utilities.regRoot;

    extends DynamicVCC.Components.Pipes.BaseClasses.FlowModels.PartialStaggeredGridFlowModel;

    replaceable package WallFriction=DynamicVCC.Components.Pipes.BaseClasses.WallFriction.Constant;

    parameter Real lambda0=1;

  equation

     if from_dp then
      m_flows=WallFriction.massFlowRate_dp(
        dps_f,
        rho_act,
        mu_act,
        lengths,
        diameters,
        (crossAreas[1:n-1]+crossAreas[2:n])/2,
        dp_small,
        lambda0);
    else
      dps_f=WallFriction.pressureLoss_m_flow(
        m_flows,
        rho_act,
        mu_act,
        lengths,
        diameters,
        (crossAreas[1:n-1]+crossAreas[2:n])/2,
        m_flow_small,
        lambda0);
     end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConstantFrictionFlow;

  model ConstantFrictionFlowPhaseChange "Constant friction factor for each phase"

    import Modelica.Fluid.Utilities.regStep;
    import Modelica.Fluid.Utilities.regRoot;
    import DynamicVCC.Media.Utilities.phaseTransition;
    extends DynamicVCC.Components.Pipes.BaseClasses.FlowModels.PartialStaggeredGridFlowModel;

    parameter Real lambda_f=1 "Liquid friction factor";

    parameter Real lambda_tp=1 "Two-phase friction factor";

    parameter Real lambda_g=1 "Vapor friction factor";

  protected
    Real lambda[n-1](each start=1);
    Medium.AbsolutePressure p[n]=Medium.pressure(states);
    Medium.AbsolutePressure p_act[n-1];
    Medium.SpecificEnthalpy h[n]=Medium.specificEnthalpy(states);
    Medium.SpecificEnthalpy h_act[n-1];
    Medium.SpecificEnthalpy hl[n-1]=Medium.bubbleEnthalpy(Medium.setSat_p(p_act));
    Medium.SpecificEnthalpy hv[n-1]=Medium.dewEnthalpy(Medium.setSat_p(p_act));
    Real x[n-1] "Thermodynamic equilibrium quality";
    Real w1[n-1], w2[n-1], w3[n-1];
  equation
    if not enableReverseFlow then
      p_act=p[1:n-1];
      h_act=h[1:n-1];
    elseif not upStreamProperties then
      p_act=0.5*(p[1:n-1]+p[2:n]);
      h_act=0.5*(h[1:n-1]+h[2:n]);
    else
      for i in 1:n-1 loop
        p_act[i]=noEvent(if m_flows[i]>0 then p[i] else p[i+1]);
        h_act[i]=noEvent(if m_flows[i]>0 then h[i] else h[i+1]);
      end for;
    end if;

    for i in 1:n-1 loop
      x[i]=(h_act[i]-hl[i])/(hv[i]-hl[i]);
      (w1[i],w2[i],w3[i])=phaseTransition(x[i]);
      lambda[i]=w1[i]*lambda_f+w2[i]*lambda_tp+w3[i]*lambda_g;
      //lambda[i]=smooth(1,noEvent(
      //if x[i]<0.5 then regStep(x[i],lambda_tp,lambda_f,1) else
      //  regStep(x[i]-1,lambda_g,lambda_tp,1)));

    end for;

    if from_dp then
      m_flows={crossAreas[i]*sqrt(2*dimensions[i]*rho_act[i]/lengths[i]/lambda[i])*regRoot(dps_f[i]) for i in 1:n-1};
    else
      dps_f={lambda[i]*lengths[i]/dimensions[i]*m_flows[i]^2/(2*rho_act[i]*crossAreas[i]^2) for i in 1:n-1};
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConstantFrictionFlowPhaseChange;

  model NominalFrictionFlow "Using nominal condition to approximate actual flow conditions"
    extends DynamicVCC.Components.Pipes.BaseClasses.FlowModels.PartialStaggeredGridFlowModel;

    replaceable package WallFriction=DynamicVCC.Components.Pipes.BaseClasses.WallFriction.Nominal;

    parameter Real k=1 "Multiplier";
    parameter Real b=1 "power constant";

  equation
    if from_dp then
      m_flows=homotopy(actual=WallFriction.massFlowRate_dp(
        dps_f,
        rho_act,
        mu_act,
        lengths,
        diameters,
        (crossAreas[1:n-1]+crossAreas[2:n])/2,
        dp_small,
        m_flow_nominal,
        dp_nominal,
        k,
        b),
        simplified=m_flow_nominal/dp_nominal*n*dps_f);
    else
      dps_f=homotopy(actual=WallFriction.pressureLoss_m_flow(
        m_flows,
        rho_act,
        mu_act,
        lengths,
        diameters,
        (crossAreas[1:n-1]+crossAreas[2:n])/2,
        m_flow_small,
        m_flow_nominal,
        dp_nominal,
        k,
        b),
        simplified=dp_nominal/n/m_flow_nominal*m_flows);
    end if;
  initial equation
    //dp_nominal=1e3*dp_small;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end NominalFrictionFlow;
end FlowModels;
