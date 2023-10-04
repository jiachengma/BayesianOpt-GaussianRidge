within DynamicVCC.Components.Pipes.BaseClasses;
partial model PartialFlowCell "1-D refrigerant flow cell (momentum balances)"

    outer DynamicVCC.Components.System system;

    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;

    parameter Integer m=1 "Number of flow cells";

    input SI.Length lengths[m] "Length of discretized cell";

    parameter Boolean enableReverseFlow=system.enableReverseFlow;

    parameter DynamicVCC.Components.Types.Dynamics momentumDynamics=system.momentumDynamics;

    // Initialization
    parameter Medium.MassFlowRate m_flow_init=system.m_flow_init;

    //Conservation quantities
    SI.Momentum I[m] "Momentum of flow cell";
    Medium.MassFlowRate m_flows[m](
    each min=if enableReverseFlow then -Modelica.Constants.inf else 0,
    each start=m_flow_init) "Mass flow rate of each flow cell (boundary of each volume cell)";

    //Source terms
    SI.Force I_flows[m] "Momentum flow across cells";
    SI.Force F_p[m] "Pressure forces";
    SI.Force F_f[m] "Frictional forces";

equation

    I={m_flows[i]*lengths[i] for i in 1:m};

    //Momentum balances
    if momentumDynamics == DynamicVCC.Components.Types.Dynamics.SteadyState then
      zeros(m) = I_flows - F_p - F_f;
    else
      lengths.*der(m_flows) = I_flows - F_p - F_f;
    end if;

initial equation
  if momentumDynamics == DynamicVCC.Components.Types.Dynamics.SteadyState_init then
    der(m_flows)=zeros(m);
  elseif momentumDynamics == DynamicVCC.Components.Types.Dynamics.Fixed_init then
    m_flows=fill(m_flow_init,m);
  end if;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end PartialFlowCell;
