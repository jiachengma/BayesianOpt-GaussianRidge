within DynamicVCC.Components.Pipes.BaseClasses;
model StaggeredGridMomentum "Staggered grid flow models (momentum grid)"

  import Modelica.Fluid.Utilities.regStep;

  // extendisng flow cells
  extends DynamicVCC.Components.Pipes.BaseClasses.PartialFlowCell(final m=n - 1);

  parameter Integer n=2 "Number of control volume";

  // Inputs to flow cells
  input Medium.ThermodynamicState states[n];

  input SI.Velocity v[n] "Flow velocity at cell boundaries";

  input SI.Area Ac[n] "Cross flow area";

  input SI.Length dimensions[n] "Characteristic dimension of fluid flow, typically pipe diameter";

  // Numerical settings
  parameter Boolean use_I_flows=false
  "= true to evaluate momentum flow differences in momentum balances";
  parameter Boolean upStreamProperties=true
  "= true to use upstream thermodynamic properties, otherwise use average properties";

  // Variables
  Medium.Density rho[n]=Medium.density(states);
  Medium.Density rho_act[n-1] "Actual density of flow cell";
  Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states);
  Medium.DynamicViscosity mu_act[n-1] "Actual viscosity of flow cell";
  Medium.AbsolutePressure p_act[n-1] "Actual pressure of flow cell";
  Medium.SaturationProperties sat[n-1]=Medium.setSat_p(p_act);

  // Friction model
  replaceable model FrictionalDP=DynamicVCC.Components.Pipes.BaseClasses.Friction.Correlations.Constant constrainedby DynamicVCC.Components.Pipes.BaseClasses.Friction.PartialFrictionalPressureDrop;

  FrictionalDP frictionalDP(
  redeclare final package Medium=Medium,
  final m=n-1,
  final rho=rho_act,
  final mu=mu_act,
  final sat=sat,
  final v=0.5*(v[1:n-1]+v[2:n]),
  final L=lengths,
  final dimension=0.5*(dimensions[1:n-1]+dimensions[2:n]));

equation

  // Properties
  if not EnableReverseFlow then
    rho_act=rho[1:n-1];
    mu_act=mu[1:n-1];
    p_act=0.5*(Medium.pressure(states[1:n-1])+Medium.pressure(states[2:n]));
  elseif not upStreamProperties then
    rho_act=0.5*(rho[1:n-1]+rho[2:n]);
    mu_act=0.5*(mu[1:n-1]+mu[2:n]);
    p_act=0.5*(Medium.pressure(states[1:n-1])+Medium.pressure(states[2:n]));
  else
    for i in 1:n-1 loop
      rho_act[i]=noEvent(if m_flows[i]>0 then rho[i] else rho[i+1]);
      mu_act[i]=noEvent(if m_flows[i]>0 then mu[i] else mu[i+1]);
      p_act[i]=noEvent(if m_flows[i]>0 then Medium.pressure(states[i]) else Medium.pressure(states[i+1]));
    end for;
  end if;

  // Source terms of momentum balances
  if use_I_flows then
    I_flows={rho[i]*v[i]^2*Ac[i] - rho[i+1]*v[i+1]^2*Ac[i+1] for i in 1:n-1};
  else
    I_flows=zeros(n-1);
  end if;

  F_p={0.5*(Ac[i]+Ac[i+1])*(Medium.pressure(states[i+1])-Medium.pressure(states[i])) for i in 1:n-1};

  //F_f={0.5*(Ac[i]+Ac[i+1])*frictionalDP.dp[i]*regStep(m_flows[i]*m_scale,1,-1) for i in 1:n-1};

  F_f=homotopy(
  actual={0.5*(Ac[i]+Ac[i+1])*frictionalDP.dp[i]*regStep(m_flows[i],1,-1) for i in 1:n-1},
  simplified=0.5*(Ac[1:n-1]+Ac[2:n])*frictionalDP.dp_nominal/frictionalDP.m_flow_nominal*m_flows);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end StaggeredGridMomentum;
