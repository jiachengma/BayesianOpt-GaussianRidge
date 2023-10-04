within DynamicVCC.Test.Chiller;
function controlSignal

  input Boolean open;
  input Real gamma_pre;
  input Real actionTime_pre;
  output Real gamma;
  output Real actionTime;

protected
  parameter Real gamma_max=1;
  parameter Real gamma_min=0.05;
  parameter Real upmax=1/150;
  parameter Real dnmax=1/38;//1/45;
  parameter Real waitTime=15;
algorithm
  if actionTime_pre<1 then
    if open and gamma_pre<gamma_max then
      gamma:=gamma_pre + upmax*actionTime_pre;
    elseif (not open) and gamma_pre>gamma_min then
      gamma:=gamma_pre - dnmax*actionTime_pre;
    else
      gamma:=gamma_pre;
    end if;
    actionTime:= 0;
  else
    if open and gamma_pre<gamma_max then
      gamma:=gamma_pre + upmax;
    elseif (not open) and gamma_pre>gamma_min then
      gamma:=gamma_pre - dnmax;
    else
      gamma:=gamma_pre;
    end if;
    actionTime:=actionTime_pre - 1;
  end if;

end controlSignal;
