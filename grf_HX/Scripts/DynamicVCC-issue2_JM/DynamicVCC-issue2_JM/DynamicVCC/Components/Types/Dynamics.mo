within DynamicVCC.Components.Types;
type Dynamics = enumeration(
    DynamicFree_init "Dynamic balance, initial guess value",
    Fixed_init "Dynamic balance, fixed inital value",
    SteadyState     "Steady-state balance",
    SteadyState_init "Dynamic balance, steady-state initialization");
