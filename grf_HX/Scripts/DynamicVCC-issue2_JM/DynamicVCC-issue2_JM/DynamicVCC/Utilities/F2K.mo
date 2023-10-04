within DynamicVCC.Utilities;
function F2K
  input SI.Temperature Tin;
  output SI.Temperature Tout;
algorithm
Tout:=(Tin - 32)*5/9 + 273.15;
end F2K;
