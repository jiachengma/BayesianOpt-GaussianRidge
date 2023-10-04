within DynamicVCC.Components.Pipes.BaseClasses;
partial model PartialVolumeCell "1-D refrigerant volume cell (mass & energy balances)"

    outer DynamicVCC.Components.System system;

    import DynamicVCC.Components.Types.DifferentialState;
    import DynamicVCC.Components.Types.Dynamics;
    import Modelica.Fluid.Utilities.regStep;
    import DynamicVCC.Media.Utilities.phaseTransition;

    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;

    parameter Integer n=2 "Number of control volumes";
    parameter Types.DifferentialState differentialState=DifferentialState.pdh;
    parameter Types.Dynamics energyDynamics=system.energyDynamics;
    parameter Types.Dynamics massDynamics=system.massDynamics;

    input SI.Volume fluidVolumes[n] "Volume of discretized cell";

/***************** Initialization ********************/
    parameter Medium.AbsolutePressure p_a_start=Medium.p_default;
    parameter Medium.AbsolutePressure p_b_start=p_a_start;
    parameter Medium.AbsolutePressure p_init[n]=linspace(p_a_start,p_b_start,n);
    parameter Medium.SpecificEnthalpy h_init[n]=fill(Medium.h_default,n);

/***************** Thermodynamic properties ********************/

    Medium.ThermodynamicState states[n];
    Medium.AbsolutePressure p[n](start=p_init);
    Medium.SpecificEnthalpy h[n](start=h_init);
    Medium.SpecificEnthalpy h_flow[n](start=h_init) "Flow-weighted enthalpy at interfaces";
    Medium.Density rho[n](start=Medium.density(Medium.setState_ph(p_init,h_init)));
    Medium.SpecificInternalEnergy u[n](start=Medium.specificInternalEnergy(Medium.setState_ph(p_init,h_init)));
    Medium.DerDensityByPressure drhodp_h[n]=Medium.density_derp_h(states)
      "partial derivative of density wrt pressure";
    Medium.DerDensityByEnthalpy drhodh_p[n]=Medium.density_derh_p(states)
      "partial derivative of density wrt enthalpy";
    Medium.MassFraction x[n]=Medium.vapourQuality(states)
    "Flow quality, equal to thermodynamic equilibrium quality in two-phase region";
    Medium.MassFraction x_s[n] "Static quality";
    Real xe[n] "Thermodynamic equilibrium quality";
    Medium.Density rho_f[n]=Medium.bubbleDensity(Medium.setSat_p(Medium.pressure(states)));
    Medium.Density rho_g[n]=Medium.dewDensity(Medium.setSat_p(Medium.pressure(states)));
    Medium.SpecificEnthalpy h_f[n]=Medium.bubbleEnthalpy(Medium.setSat_p(Medium.pressure(states)));
    Medium.SpecificEnthalpy h_g[n]=Medium.dewEnthalpy(Medium.setSat_p(Medium.pressure(states)));

    // Conservation quantities
    SI.MassFlowRate dMdt[n];
    SI.EnergyFlowRate dUdt[n];

    // Source terms
    Medium.MassFlowRate m_flows_cell[n] "Mass flow rate across each cell";
    SI.EnthalpyFlowRate H_flows_cell[n] "Enthalpy flow rate across each cell";
    SI.HeatFlowRate Q_flows[n] "Heat flow rate of each cell";

    // Void fraction
    Real gamma[n](each min=0, each max=1) "Void fraction";
    replaceable model SlipRatio=DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.Homogeneous
     constrainedby DynamicVCC.Components.Pipes.BaseClasses.SlipRatio.PartialSlipRatio;

    SlipRatio slipRatio(
      redeclare package Medium = Medium,
      final n=n,
      final states=states);

protected
    SI.DerEnthalpyByPressure dhdp_rho[n] "partial derivative of enthalpy wrt pressure at constant density";
    SI.DerDensityByPressure drhodp_u[n] "partial derivative of density wrt pressure at constant internal energy";
    SI.DerEnthalpyByPressure dudp_rho[n] "partial derivative of internal energy wrt pressure at constant density";
    SI.SpecificEnthalpy h_fg[n]=h_g-h_f;
    Real w_l[n](each min=0,each max=1) "liquid-phase weights";
    Real w_tp[n](each min=0,each max=1) "two-phase weights";
    Real w_v[n](each min=0,each max=1) "vapor-phase weights";
equation

    // Select state variables
    if differentialState==DifferentialState.ph then
      dMdt={fluidVolumes[i]*(drhodp_h[i]*der(p[i]) + drhodh_p[i]*der(h[i])) for i in 1:n};
      dUdt={fluidVolumes[i]*(h[i]*drhodp_h[i] - 1)*der(p[i]) + fluidVolumes[i]*(h[i]*drhodh_p[i] + rho[i])*der(h[i]) for i in 1:n};
      rho={(w_l[i]+w_v[i])*Medium.density(states[i])+w_tp[i]*(rho_g[i]*gamma[i]+rho_f[i]*(1-gamma[i])) for i in 1:n};
      u=Medium.specificInternalEnergy(states);
      gamma={x[i]/(x[i]+rho_g[i]/rho_f[i]*slipRatio.S[i]*(1-x[i])) for i in 1:n};
    elseif differentialState==DifferentialState.pdh then
      dMdt={fluidVolumes[i]*der(rho[i]) for i in 1:n};
      dUdt={fluidVolumes[i]*((rho[i]*dhdp_rho[i]-1)*der(p[i])+(rho[i]/drhodh_p[i]+h[i])*der(rho[i])) for i in 1:n};
      der(h)={dhdp_rho[i]*der(p[i])+der(rho[i])/drhodh_p[i] for i in 1:n};
      //h={w_l[i]*Medium.liquidSpecificEnthalpy_pd(p[i],rho[i])+w_v[i]*Medium.vaporSpecificEnthalpy_pd(p[i],rho[i])+w_tp[i]*(x_s[i]*h_g[i]+(1-x_s[i])*h_f[i]) for i in 1:n};
      u=Medium.specificInternalEnergy(states);
      gamma={w_tp[i]*(rho[i]-rho_f[i])/(rho_g[i]-rho_f[i])+w_v[i]*1 for i in 1:n};
    elseif differentialState==DifferentialState.pdu then
      dMdt={fluidVolumes[i]*der(rho[i]) for i in 1:n};
      dUdt={fluidVolumes[i]*(rho[i]*der(u[i])+u[i]*der(rho[i])) for i in 1:n};
      der(p)={der(rho[i])/drhodp_u[i]+der(u[i])/dudp_rho[i] for i in 1:n};
      h={u[i]+p[i]/rho[i] for i in 1:n};
      gamma={w_tp[i]*(rho[i]-rho_f[i])/(rho_g[i]-rho_f[i])+w_v[i]*1 for i in 1:n};
/*
    elseif differentialState==DifferentialState.du then
      dMdt={V[i]*der(rho[i]) for i in 1:n};
      dUdt={V[i]*(rho[i]*der(u[i])+u[i]*der(rho[i])) for i in 1:n};
      p={Medium.pressure_du(rho[i],u[i]) for i in 1:n};
      h={u[i]+p[i]/rho[i] for i in 1:n};
*/
    else
      assert(false, "Unknown differential states");
    end if;

    // Use p,h inputs for properties evaluation
    xe=(h-h_f)./(h_g-h_f);
    for i in 1:n loop
      states[i] = Medium.setState_ph(p[i],h[i]);
      x_s[i]=x[i]/(x[i]+slipRatio.S[i]*(1-x[i]));
      (w_l[i],w_tp[i],w_v[i])=phaseTransition(xe[i]);
      h_flow[i]=h[i]+(x[i]-x_s[i])*h_fg[i];
      dhdp_rho[i]=-drhodp_h[i]/drhodh_p[i];
      drhodp_u[i]=(drhodp_h[i]+drhodh_p[i]/rho[i])/(1+p[i]/rho[i]^2*drhodh_p[i]);
      dudp_rho[i]=-(drhodp_h[i]+drhodh_p[i]/rho[i])/drhodh_p[i];
    end for;

    // Mass balances
    dMdt = m_flows_cell;

    // Energy balances
    dUdt = H_flows_cell + Q_flows;

initial equation

  if energyDynamics==Dynamics.SteadyState_init then
    if differentialState==DifferentialState.ph then
      der(h)=zeros(n);
    elseif differentialState==DifferentialState.pdh then
      der(h)=zeros(n);
    elseif differentialState==DifferentialState.pdu then
      der(u)=zeros(n);
    end if;
  elseif energyDynamics==Dynamics.Fixed_init then
    if differentialState==DifferentialState.ph then
      h=h_init;
    elseif differentialState==DifferentialState.pdh then
      h=h_init;
      //rho=Medium.density(Medium.setState_ph(p_init,h_init));
    elseif differentialState==DifferentialState.pdu then
      u=Medium.specificInternalEnergy(Medium.setState_phX(p_init,h_init));
    end if;
  end if;

  if massDynamics==Dynamics.SteadyState_init then
    der(p)=zeros(n);
  elseif massDynamics==Dynamics.Fixed_init then
    p=p_init;
  end if;


  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end PartialVolumeCell;
