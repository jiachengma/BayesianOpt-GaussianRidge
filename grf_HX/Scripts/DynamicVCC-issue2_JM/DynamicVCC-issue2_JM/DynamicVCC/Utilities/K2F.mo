within DynamicVCC.Utilities;
function K2F
  input SI.Temperature Tin;
  output SI.Temperature Tout;
algorithm
Tout:=(Tin - 273.15)*9/5 + 32;
end K2F;
