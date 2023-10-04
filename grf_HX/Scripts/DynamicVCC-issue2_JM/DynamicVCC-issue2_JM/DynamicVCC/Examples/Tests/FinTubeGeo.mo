within DynamicVCC.Examples.Tests;
record FinTubeGeo "Calculate geometry parameters of fin-tube heat exchanger."
  extends Modelica.Icons.Record;

  import Modelica.Constants.pi;
  parameter SI.Diameter D_o=1 "tube outter diameter";
  parameter SI.Diameter D_i=1;
  parameter SI.SpecificHeatCapacity cp_tube=385;
  parameter SI.Density rho_tube=8900;
  parameter SI.SpecificHeatCapacity cp_fin=900;
  parameter SI.Density rho_fin=2700;
  parameter SI.Length L_tube=1 "Length of single tube";
  parameter Integer N_tuberow=1 "Number of tube rows";
  parameter Integer N_t_prow=1 "Number of tubes per row";
  parameter Integer N_circuits=1 "Number of circuits";
  parameter SI.Length pf=1 "distance between fins";
  parameter SI.Length pt=1 "Tube transverse spacing";
  parameter SI.Length pl=1 "Tube longitudinal spacing";
  parameter Real fin_meter=1 "Number of fins per meter";
  parameter SI.Thickness t_fin=1 "Fin thickness";

  parameter SI.Area Ac_tube=pi*D_i^2/4;
  parameter SI.Area Ao_tube=pi*D_o*L_tube "tube outter surface area";
  parameter SI.Area Ai_tube=pi*D_i*L_tube "tube inner surface area";
  parameter SI.Volume V_tube=Ac_tube*L_tube;
  parameter Integer N_tube_tot=N_tuberow*N_t_prow;
  parameter SI.Mass M_wall=pi/4*(D_o^2-D_i^2)*L_tube*rho_tube "Mass of each tube";

  parameter SI.Area HTA_r=Ai_tube*N_tube_tot;
  parameter SI.Area Ac_r=Ac_tube*N_circuits;
  parameter SI.Volume V_r=V_tube*N_tube_tot;
  parameter SI.Mass M_tube=M_wall*N_tube_tot;
  parameter SI.Length L_circuit=N_tube_tot*L_tube/N_circuits "Average circuit length";
  parameter Real secant=sqrt(pl^2+(pt/2)^2)/pl;

  parameter Real Eta_fin_overall=1;
  parameter SI.Area Atube=N_tuberow*N_t_prow*pi*D_o*L_tube "Total tube outter surface area";
  parameter SI.Height H_hx=N_t_prow*pt "Heat exchanger height";
  parameter Real N_fin_ptube=L_tube*fin_meter "Number of fins along the tube";
  parameter SI.Area Ac_e=H_hx*L_tube-t_fin*N_fin_ptube*(H_hx-D_o*N_t_prow)-N_t_prow*D_o*L_tube "Coil cross sectional area";
  parameter SI.Area HTA_fin=N_fin_ptube*2*(H_hx*pl*N_tuberow*secant-N_tuberow*N_t_prow*pi*D_o^2/4);
  parameter SI.Area HTA_tube=N_tuberow*N_t_prow*pi*D_o*(L_tube-L_tube*fin_meter*t_fin);
  parameter SI.Area HTA_e=HTA_fin+HTA_tube "Total air side heat transfer area";
  parameter SI.Mass M_fin=N_fin_ptube*(H_hx*pl*N_tuberow-N_tuberow*N_t_prow*pi*D_o^2/4)*t_fin*rho_fin "Total fin mass";

end FinTubeGeo;
