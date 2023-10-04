within DynamicVCC.Components.Types;
type DifferentialState = enumeration(
    ph
   "pressure and enthalpy",
    pdh
    "pressure, density and enthalpy",
    pdu
    "pressure, density and internal energy",
    du
    "density and internal energy")
  "Select state variables and differential formula for conservation terms";
