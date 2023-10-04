within DynamicVCC.Test.Chiller;
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
