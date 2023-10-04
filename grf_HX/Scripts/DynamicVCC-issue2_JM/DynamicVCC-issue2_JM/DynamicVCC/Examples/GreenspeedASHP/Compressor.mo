within DynamicVCC.Examples.GreenspeedASHP;
model Compressor

  import Modelica.Fluid.Utilities.regPow;
  extends DynamicVCC.Components.Units.MassFlowDevices.Compressor.PartialEfficiencyMap;

  parameter SI.ThermalConductance UA=1 "Shell wall thermal conductance";
  parameter Real n=1.3 "Polytropic coefficient";

protected
  Real eta_v(min=0,max=1,start=0.8);
  Real eta_p(min=0,max=1,start=0.7);
  SI.Temperature T_w "Average shell wall temperature";
equation

  // Efficiency maps
  eta_v=-0.00146604*speed-0.0029233*p_ratio+0.92979595;

  eta_p=0.08615504*p_ratio+0.01499364*speed-0.00288711*p_ratio*speed+0.15308804;

  // Flow rate
  m_flow=homotopy(actual=eta_v*Vs*rho_suc*speed,
                simplified=m_flow_nominal*speed/speed_nominal);
  // Power
  Pwr=m_flow*(n/(n - 1)*p_suc/rho_suc*(regPow(p_ratio,(n-1)/n) - 1))/eta_p;

  // Heat loss
  T_w=T_dis;
  Q_loss=UA*(T_w-T_amb);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
end Compressor;
