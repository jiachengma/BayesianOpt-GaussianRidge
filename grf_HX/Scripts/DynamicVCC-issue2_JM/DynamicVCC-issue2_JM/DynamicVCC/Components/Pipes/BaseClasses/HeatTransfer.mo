within DynamicVCC.Components.Pipes.BaseClasses;
package HeatTransfer
  extends Modelica.Icons.VariantsPackage;
  partial model PartialFlowHeatTransfer
    extends DynamicVCC.Interfaces.PartialHeatTransfer;

    input SI.Length dimensions[n] "Characteristic dimension, typically plate length or tube diameter";

    input SI.Length lengths[n] "Length along flow path";

    input SI.Velocity vs[n] "Mean velocity of each control volume";

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
            extent={{-80,80},{80,-80}},
            lineColor={238,46,47},
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid)}),                      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialFlowHeatTransfer;

  partial model PartialHeatTransferCoefficient "Partial model for calculating heat transfer coefficient, no heatport"

    replaceable package Medium=Modelica.Media.Interfaces.PartialTwoPhaseMedium;

    parameter Integer n=1 "Number of control volumes";

    input Medium.ThermodynamicState states[n] "Fluid states";

    input SI.Area surfaceAreas[n] "Heat transfer area";

    input SI.Length dimensions[n] "Characteristic dimension, typically plate length or tube diameter";

    input SI.Length lengths[n] "Length along flow path";

    input SI.Velocity vs[n] "Mean velocity of each control volume";

    parameter SI.CoefficientOfHeatTransfer alpha0=100 "Guess values for heat transfer coefficients";

    output SI.CoefficientOfHeatTransfer alphas[n](each start=alpha0) "Heat transfer coefficient";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
            extent={{-80,80},{80,-80}},
            lineColor={238,46,47},
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid)}),                      Diagram(coordinateSystem(preserveAspectRatio=false)));
  end PartialHeatTransferCoefficient;

  model ConstantFlowHeatTransfer "Constant heat transfer coefficient for calculating heat flow"
    extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer;

    parameter SI.CoefficientOfHeatTransfer alpha0=1 "heat transfer coefficient";

  equation
    Q_flows={alpha0*surfaceAreas[i]*(heatPorts[i].T-T[i]) for i in 1:n};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConstantFlowHeatTransfer;

  model ConstantFlowPhaseChange "Constant heat transfer coefficient for each phase"

    import Modelica.Fluid.Utilities.regStep;

    extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer;

    parameter SI.CoefficientOfHeatTransfer alpha_f=100 "Liquid HTC";

    parameter SI.CoefficientOfHeatTransfer alpha_tp=100 "Two-phase HTC";

    parameter SI.CoefficientOfHeatTransfer alpha_g=100 "Vapor HTC";

  protected
    SI.CoefficientOfHeatTransfer alphas[n];
    Medium.SaturationProperties sat[n]=Medium.setSat_p(Medium.pressure(states));
    Medium.SpecificEnthalpy hl[n]=Medium.bubbleEnthalpy(sat);
    Medium.SpecificEnthalpy hv[n]=Medium.dewEnthalpy(sat);
    Real x[n] "Thermodynamic equilibrium quality";

  equation
    x={(Medium.specificEnthalpy(states[i])-hl[i])/(hv[i]-hl[i]) for i in 1:n};
    alphas={smooth(1,noEvent(if x[i]<0.5 then regStep(x[i],alpha_tp,alpha_f,1e-4) else
    regStep(x[i]-1,alpha_g,alpha_tp,1e-4))) for i in 1:n};
    Q_flows={alphas[i]*surfaceAreas[i]*(heatPorts[i].T-T[i]) for i in 1:n};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConstantFlowPhaseChange;

  model ConstantHeatTransfer "Constant heat transfer coefficient"
    extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialHeatTransferCoefficient;

  equation
    alphas=fill(alpha0,n);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ConstantHeatTransfer;

  package Correlations
    extends Modelica.Icons.VariantsPackage;
    partial model PartialPipeHeatTransferCoefficient "Base class for pipe heat transfer correlation"
      import Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers;

      input SI.Temperature Ts[n] "Wall surface temperatures";

      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialHeatTransferCoefficient;

      SI.ReynoldsNumber Re[n];
      SI.PrandtlNumber Pr[n]=Medium.prandtlNumber(states);
      SI.NusseltNumber Nu[n];
      Medium.Density d[n]=Medium.density(states);
      Medium.Temperature T[n]=Medium.temperature(states);
      Medium.ThermalConductivity lambda[n]=Medium.thermalConductivity(states);
      Medium.DynamicViscosity mu[n]=Medium.dynamicViscosity(states);
    equation
      Re=CharacteristicNumbers.ReynoldsNumber(vs, d, mu, dimensions);
      Nu=CharacteristicNumbers.NusseltNumber(alphas, dimensions, lambda);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end PartialPipeHeatTransferCoefficient;

    model SinglePhase_Gnielinski "Gnielinski equation for single phase flow (VDI Heat Atlas 2010)"

      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.PartialPipeHeatTransferCoefficient;

    protected
      parameter SI.ReynoldsNumber Re_lam=2300;
      parameter SI.ReynoldsNumber Re_tur=10000;
      Real Num_1;
      Real Num_2[n];
      Real Num_3[n];
      Real Nu_lam[n] "Mean Nu for laminar flow Re<2300";
      Real Nu_turb[n] "Mean Nu for turbulent flow Re>10000";
      Real Nu_trans[n] "Mean Nu in transition region";
      Real f[n] "Friction factor";
      Real gamma[n]
      "Interpolation for the transition regime between laminar and turbulent";

    equation
      Num_1=3.66;
      for i in 1:n loop
        Num_2[i]=1.615*(Re[i]*Pr[i]*dimensions[i]/lengths[i])^(1/3);
        Num_3[i]=(2/(1+22*Pr[i]))^(1/6)*(Re[i]*Pr[i]*dimensions[i]/lengths[i])^(1/2);
        Nu_lam[i]=(Num_1^3+0.7^3+(Num_2[i]-0.7)^3+Num_3[i]^3)^(1/3);
        f[i]=(1.8*log10(Re[i])-1.5)^(-2);
        Nu_turb[i]=(f[i]/8*Re[i]*Pr[i])/(1+12.7*sqrt(f[i]/8)*(Pr[i]^(2/3)-1))*(1+(dimensions[i]/lengths[i])^(2/3));
        gamma[i]=(Re[i]-Re_lam)/(Re_tur-Re_lam);
        Nu_trans[i]=(1-gamma[i])*Nu_lam[i]+gamma[i]*Nu_turb[i];
        Nu[i]=smooth(0,noEvent(if Re[i]<Re_lam then Nu_lam[i] elseif Re[i]<Re_tur then Nu_trans[i] else Nu_turb[i]));
      end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end SinglePhase_Gnielinski;

    model Boiling_GungorWinterton "Gungor and Winterton (1987)"
      import Modelica.Fluid.Utilities.regPow;
      import Modelica.Fluid.Utilities.regStep;
     extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.PartialPipeHeatTransferCoefficient;

     parameter SI.Time tau=1e-3 "Time constant";
    protected
      parameter Real eps=1e-6;
      Medium.SaturationProperties sat[n]=Medium.setSat_p(Medium.pressure(states));
      Medium.ThermodynamicState bubbleStates[n];
      Medium.ThermodynamicState dewStates[n];
      Medium.Density d_l[n]=Medium.density(bubbleStates);
      Medium.Density d_v[n]=Medium.density(dewStates);
      Medium.DynamicViscosity mu_l[n]=Medium.dynamicViscosity(bubbleStates);
      Medium.DynamicViscosity mu_v[n]=Medium.dynamicViscosity(dewStates);
      Medium.PrandtlNumber Pr_l[n]=Medium.prandtlNumber(bubbleStates);
      Medium.ThermalConductivity k_l[n]=Medium.thermalConductivity(bubbleStates);
      SI.SpecificEnthalpy h_fg[n]=Medium.specificEnthalpy(dewStates)-Medium.specificEnthalpy(bubbleStates);
      Medium.MassFraction x[n]=Medium.vapourQuality(states);
      SI.ReynoldsNumber Re_l[n] "liquid only flow Reynolds Number";
      SI.CoefficientOfHeatTransfer alphas_l[n] "liquid only heat transfer coefficient";
      Real Bo[n] "Boiling number";
      SI.FroudeNumber Fr[n];
      SI.HeatFlux q[n];
      Real E[n] "enhancement factor";
      Real S[n] "suppression factor";
      Real E2[n] "multiplication for horizontal tube and Fr<0.05";
      Real S2[n] "multiplication for horizontal tube and Fr<0.05";
      SI.CoefficientOfHeatTransfer alphas_nominal[n](each start=alpha0);

    equation
      Re_l=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(vs,d_l,mu_l,dimensions);
      for i in 1:n loop
        q[i]=alphas[i]*(Ts[i]-T[i]);
        bubbleStates[i]=Medium.setBubbleState(sat[i]);
        dewStates[i]=Medium.setDewState(sat[i]);
        Bo[i]=noEvent(abs(q[i]))/(h_fg[i]*d[i]*vs[i]);
        Fr[i]=(d[i]*vs[i]+eps)^2/(d_l[i]^2*9.81*dimensions[i]);
        //E2[i]=if Fr[i]<0.05 then Fr[i]^(0.1-2*Fr[i]) else 1;
        //S2[i]=if Fr[i]<0.05 then sqrt(Fr[i]) else 1;
        E2[i]=regStep(Fr[i]-0.05,1,Fr[i]^(0.1-2*Fr[i]),1e-6);
        S2[i]=regStep(Fr[i]-0.05,1,sqrt(Fr[i]),1e-6);
        E[i]=1.12*(x[i]/(1-x[i]+eps))^0.75*(d_l[i]/d_v[i])^0.41;
        S[i]=1+3e3*regPow(Bo[i],0.86,1e-5);
        alphas_l[i]=0.023*Re_l[i]^0.8*Pr_l[i]^0.4*k_l[i]/dimensions[i];
        alphas_nominal[i]=alphas_l[i]*(S[i]*S2[i]+E[i]*E2[i]);
      end for;

      der(alphas)=(alphas_nominal-alphas)/tau;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Boiling_GungorWinterton;

    model Condensation_Shah "Film condensation correlation (Shah, 1979)"
      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.PartialPipeHeatTransferCoefficient;

      parameter SI.AbsolutePressure pc=Medium.fluidConstants[1].criticalPressure "Critical pressure";

      parameter SI.Time tau=1e-3 "Time constant";

    protected
      Medium.SaturationProperties sat[n]=Medium.setSat_p(Medium.pressure(states));
      Medium.ThermodynamicState bubbleStates[n];
      Medium.Density d_l[n]=Medium.density(bubbleStates);
      Medium.DynamicViscosity mu_l[n]=Medium.dynamicViscosity(bubbleStates);
      Medium.PrandtlNumber Pr_l[n]=Medium.prandtlNumber(bubbleStates);
      Medium.ThermalConductivity k_l[n]=Medium.thermalConductivity(bubbleStates);
      Medium.MassFraction x[n]=Medium.vapourQuality(states);
      SI.ReynoldsNumber Re_l[n] "Assume all the mass flowing as liquid";
      SI.CoefficientOfHeatTransfer alphas_l[n] "Assume all the mass flowing as liquid";
      Real p_r[n] "Reduced pressure=p/pc";
      SI.CoefficientOfHeatTransfer alphas_nominal[n](each start=alpha0);
    equation
      Re_l=Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(vs,d_l,mu_l,dimensions);
      p_r=Medium.pressure(states)/pc;
      for i in 1:n loop
        bubbleStates[i]=Medium.setBubbleState(sat[i]);
        alphas_l[i]=0.023*Re_l[i]^0.8*Pr_l[i]^0.4*k_l[i]/dimensions[i];
        alphas_nominal[i]=alphas_l[i]*((1-x[i])^0.8+3.8*x[i]^0.76*(1-x[i])^0.04/p_r[i]^0.38);
      end for;

      der(alphas)=(alphas_nominal-alphas)/tau;


      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Condensation_Shah;

    model NominalHeatTransfer "Using nominal condition heat transfer coefficient"
      import DynamicVCC.Utilities.regPowGen;
      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.PartialPipeHeatTransferCoefficient;

      parameter SI.CoefficientOfHeatTransfer alpha_nominal=1e3;
      parameter SI.MassFlowRate m_flow_nominal=1;
      parameter Real k=1 "Multiplier";
      parameter Real b=1 "power constant";
      parameter SI.Area crossAreas[n];

    protected
      Medium.MassFlowRate m_flows[n];
      Medium.Density rhos[n]=Medium.density(states);
    equation
      m_flows={rhos[i]*vs[i]*crossAreas[i] for i in 1:n};
      alphas={k*alpha_nominal*abs(regPowGen(m_flows[i]/m_flow_nominal,b,1e-3)) for i in 1:n};

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end NominalHeatTransfer;

    model CorrelationPhaseChange "single-phase and two-phase correlations transition based on phase, Fuzzy modeling for smooth transitions"

      import Modelica.Fluid.Utilities.regStep;

      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialFlowHeatTransfer;

      replaceable model VaporPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.SinglePhase_Gnielinski;

      replaceable model LiquidPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.SinglePhase_Gnielinski;

      replaceable model TwoPhase=DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.Correlations.Condensation_Shah;

      VaporPhase vaporPhase(
      redeclare final package Medium=Medium,
      final n=n,
      final states=states,
      final surfaceAreas=surfaceAreas,
      final dimensions=dimensions,
      final lengths=lengths,
      final vs=vs,
      final Ts=heatPorts.T,
      final alpha0=alpha0);

      LiquidPhase liquidPhase(
      redeclare final package Medium=Medium,
      final n=n,
      final states=states,
      final surfaceAreas=surfaceAreas,
      final dimensions=dimensions,
      final lengths=lengths,
      final vs=vs,
      final Ts=heatPorts.T,
      final alpha0=alpha0);

      TwoPhase twoPhase(
      redeclare final package Medium=Medium,
      final n=n,
      final states=states,
      final surfaceAreas=surfaceAreas,
      final dimensions=dimensions,
      final lengths=lengths,
      final vs=vs,
      final Ts=heatPorts.T,
      final alpha0=alpha0);

      parameter SI.CoefficientOfHeatTransfer alpha0=100;
      SI.CoefficientOfHeatTransfer alphas[n](each start=alpha0);

    protected
      Medium.SpecificEnthalpy hl[n] "Bubble line enthalpy";
      Medium.SpecificEnthalpy hv[n] "Dew line enthalpy";
      Real x[n] "Thermodynamic equilibrium quality";
      Real NL[n] "Fuzzy number for liquid region";
      Real NT[n] "Fuzzy number for two-phase region";
      Real NV[n] "Fuzzy number for vapor region";
      Real w1[n] "Liquid HTC weight";
      Real w2[n] "Two-phase HTC weight";
      Real w3[n] "Vapor HTC weight";
    equation
      hl=Medium.bubbleEnthalpy(Medium.setSat_p(Medium.pressure(states)));
      hv=Medium.dewEnthalpy(Medium.setSat_p(Medium.pressure(states)));
      x=(Medium.specificEnthalpy(states)-hl)./(hv-hl);
      //alphas={smooth(1,noEvent(if x[i]<0.5 then regStep(x[i],twoPhase.alphas[i],liquidPhase.alphas[i]) else
      //regStep(x[i]-1,vaporPhase.alphas[i],twoPhase.alphas[i]))) for i in 1:n};

      // Fuzzy modeling
      for i in 1:n loop
        NL[i]=regStep(x[i],0,1,1e-4);
        NV[i]=regStep(x[i]-1,1,0,1e-4);
        NT[i]=1-NL[i]-NV[i];
      end for;
      w1=NL./(NL+NV+NT);
      w2=NT./(NL+NV+NT);
      w3=NV./(NL+NV+NT);
      alphas={w1[i]*liquidPhase.alphas[i]+w2[i]*twoPhase.alphas[i]+w3[i]*vaporPhase.alphas[i] for i in 1:n};

      Q_flows={alphas[i]*surfaceAreas[i]*(heatPorts[i].T-T[i]) for i in 1:n};


      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end CorrelationPhaseChange;

    model PlainFinAndTube_Wang "Air flow heat transfer coefficient for fin-and-tube coil (Wang, 2000)"
      // Note: Forced convection only
      import Modelica.Fluid.Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber;
      import Modelica.Fluid.Utilities.regPow;
      extends DynamicVCC.Components.Pipes.BaseClasses.HeatTransfer.PartialHeatTransferCoefficient;

      parameter SI.Time tau=1e-3 "Time constant";

      // Coil geometry
      parameter Integer nRow=1 "Number of tube rows";

      parameter SI.Length pf=1 "Fin pitch";

      parameter SI.Length pl=1 "Tube longitudinal tube pitch";

      parameter SI.Length pt=1 "Tube transverse tube pitch";

      parameter SI.Thickness t_fin=1 "Fin thickness";

      parameter SI.Diameter Dh "hydraulic diamters";

    protected
      SI.Diameter Dc[n] "fin collar outside diameters";
      Medium.Density rhos[n]=Medium.density(states);
      Medium.DynamicViscosity mus[n]=Medium.dynamicViscosity(states);
      Medium.PrandtlNumber Pr[n]=Medium.prandtlNumber(states);
      Medium.ThermalConductivity k[n]=Medium.thermalConductivity(states);
      SI.ReynoldsNumber Re_Dc[n] "Reynolds number based on fin collar outside diameters";
      SI.NusseltNumber Nu[n];
      SI.CoefficientOfHeatTransfer alphas_nominal[n](each start=alpha0);
      Real P1[n],P2[n],P3[n],P4[n],P5[n],P6[n] "Correlation parameters";
      Real j[n] "Colburn j factors";
    equation
      Dc=dimensions+2*fill(t_fin,n);
      for i in 1:n loop
        Re_Dc[i]=noEvent(max(1e-6,ReynoldsNumber(vs[i],rhos[i],mus[i],Dc[i])));
        P1[i]=1.9-0.23*log(Re_Dc[i]);
        P2[i]=-0.236+0.126*log(Re_Dc[i]);
        P3[i]=-0.361-0.042*nRow/log(Re_Dc[i])+0.158*log(nRow*(pf/Dc[i])^0.41);
        P4[i]=-1.224-0.076*(pl/Dh)^1.42/log(Re_Dc[i]);
        P5[i]=-0.083+0.058*nRow/log(Re_Dc[i]);
        P6[i]=-5.735+1.21*log(Re_Dc[i]/nRow);
        if nRow==1 then
          j[i]=0.108*regPow(Re_Dc[i],-0.29,1)*(pt/pl)^P1[i]*(pf/Dc[i])^(-1.084)*(pf/Dh)^(-0.786)*(pf/pt)^P2[i];
        else
          j[i]=0.086*regPow(Re_Dc[i],P3[i],1)*nRow^P4[i]*(pf/Dc[i])^P5[i]*(pf/Dh)^P6[i]*(pf/pt)^(-0.93);
        end if;
        Nu[i]=alphas_nominal[i]*dimensions[i]/k[i];
        j[i]=Nu[i]/(Re_Dc[i]*Pr[i]^(1/3));
      end for;

      der(alphas)=(alphas_nominal-alphas)/tau;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p><span style=\"font-size: 14pt;\">Heat transfer and friction characteristics of plain fin-and-tube heat exchangers, part II: Correlation (Wang, 2000)</span></p>
</html>"));
    end PlainFinAndTube_Wang;
  end Correlations;
end HeatTransfer;
