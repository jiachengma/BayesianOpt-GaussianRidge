within DynamicVCC.Components.Units.HX;
package MB "Moving boundary heat exchanger models"

  package Components
    package BaseClasses
      extends Modelica.Icons.BasesPackage;

      partial model RefControlVolume

        replaceable package Medium = ExternalMedia.Media.CoolPropMedium constrainedby Modelica.Media.Interfaces.PartialMedium;
          Medium.ThermodynamicState state;
          Medium.SaturationProperties satstate;
          Medium.ThermodynamicState satstate_g "Dew point";
          Medium.ThermodynamicState satstate_f "Bubble point";

        /*************** Connectors **********************/
        TransientVCC.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium);
        TransientVCC.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium);
        DynamicVCC.Components.Units.HX.MB.Interfaces.MB_HeatPort_out heatport;
        DynamicVCC.Components.Units.HX.MB.Interfaces.MB_port mb_port_a;
        DynamicVCC.Components.Units.HX.MB.Interfaces.MB_port mb_port_b;

        /*************** Geometry ***********************/
        parameter SI.Length L=1; //Total HX length
        parameter SI.Area Ac=1; //Cross-sectional area
        parameter SI.Area HTA=1 "heat transfer area";

        /************** Variables ************************/
        parameter Medium.AbsolutePressure p_scale=1;
        parameter Medium.SpecificEnthalpy h_scale=1;
        //parameter SI.CoefficientOfHeatTransfer HTC_nominal=1;
        parameter Boolean Const_HTC=true;
         parameter SI.CoefficientOfHeatTransfer HTC_nominal=1;
        parameter Medium.MassFlowRate m_nominal=1;
        parameter Medium.AbsolutePressure p_ini=Medium.reference_p;
        parameter Medium.SpecificEnthalpy h_ini=Medium.h_default;
        parameter Boolean Steadystate_ini=true;
        parameter Real zeta_a_ini=0;
        parameter Real zeta_b_ini=0;

        Real zeta_a(min=0,max=1,start=zeta_a_ini) "Normalized length at inlet";
        Real zeta_b(min=0,max=1,start=zeta_b_ini) "Normalized length at outlet";
        Real p(start=p_ini/p_scale,fixed=true); //Scaled!
        Medium.SpecificEnthalpy h;
        Medium.SpecificEnthalpy h_in "Inlet enthalpy";
        Real h_out(start=h_ini/h_scale) "Outlet enthalpy";
        Medium.Temperature T(displayUnit="degC");
        Medium.Temperature Tt(displayUnit="degC") "Tube wall temperature";
        Medium.MassFlowRate m_dot_in "Inlet mass flow rate";
        Medium.MassFlowRate m_dot_out "Outlet mass flow rate";
        SI.EnergyFlowRate Q_dot "Heat transfer rate from working fluid";
        parameter Boolean PhaseChange=true; //if true this CV is connected with others, otherwise it models the whole HX.
        parameter Boolean subcooled=true; //true if subcooled, otherwise superheated for single phase zone. For two-phase zone it indicates the single phase zone that is connected to it.
        parameter Boolean InletZone=true; //true if at the inlet of HX, false if at the outlet.

        /*
  parameter Boolean InletZone=true; //true if at the inlet of HX, false if at the outlet.
  //An exception rule is that when PhaseChange is false and InletZone is true, then subcooled=true is three-zone condenser, subcooled=false is three-zone evaporator, otherwise it is two-phase zone alone.
  */

      protected
        SI.Length L_ab "Zone length";
        Medium.DerDensityByPressure drhodp "partial derivative of density wrt pressure";
        Medium.DerDensityByEnthalpy drhodh "partial derivative of density wrt enthalpy";
        SI.CoefficientOfHeatTransfer HTC;
        Medium.Density rho_g "Dew density";
        Medium.Density rho_f "Bubble density";
        Medium.SpecificEnthalpy h_g "Dew enthalpy";
        Medium.SpecificEnthalpy h_f "Bubble enthalpy";
        Medium.DerEnthalpyByPressure dhdp_g "Partial der of dew enthalpy wrt pressure";
        Medium.DerEnthalpyByPressure dhdp_f "Partial der of bubble enthalpy wrt pressure";

      equation

        L_ab=L*(zeta_b-zeta_a);
        Q_dot = HTC_nominal * L_ab/L*HTA * (T-Tt);

        if Const_HTC then
          HTC=HTC_nominal;
        else
          HTC=HTC_nominal*((m_dot_in+m_dot_out)/2/m_nominal)^2;
        end if;

        /******************** Thermodynamic Properties *********************/
        state=Medium.setState_phX(p*p_scale,h);
        satstate=Medium.setSat_p(p*p_scale);
        satstate_g=Medium.setDewState(satstate);
        satstate_f=Medium.setBubbleState(satstate);
        T=Medium.temperature(state);
        drhodp=Medium.density_derp_h(state);
        drhodh=Medium.density_derh_p(state);
        h=(h_in+h_out*h_scale)/2;

        rho_g=satstate_g.d;
        rho_f=satstate_f.d;
        h_g=satstate_g.h;
        h_f=satstate_f.h;
        dhdp_g=Medium.dDewEnthalpy_dPressure(satstate);
        dhdp_f=Medium.dBubbleEnthalpy_dPressure(satstate);

        /***************** Boundary conditions ***********************/
        m_dot_in=inlet.m_flow;
        m_dot_out=-outlet.m_flow;
        Tt=heatport.T;
        heatport.Q_flow=-Q_dot;
        heatport.L_a=zeta_a*L;
        heatport.L_b=zeta_b*L;
        h_in=inlet.h_outflow;
        inlet.h_outflow=inStream(inlet.h_outflow);
        h_out*h_scale=outlet.h_outflow;
        inlet.p=p*p_scale;
        outlet.p=p*p_scale;
        mb_port_a.L=zeta_a*L;
        mb_port_b.L=zeta_b*L;

      initial equation
        if Steadystate_ini then
          //der(p)=0;
          der(zeta_a)=0;
          der(zeta_b)=0;
        end if;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end RefControlVolume;

      partial model MetalWall "1-D lumped capacity model"

        /*********************Connectors*********************************/
      DynamicVCC.Components.Units.HX.MB.Interfaces.MB_HeatPort_in heatport_in;
      DynamicVCC.Components.Units.HX.MB.Interfaces.MB_HeatPort_out heatport_out;

        parameter SI.Mass M=1;
        parameter SI.Length L=1;
        parameter SI.SpecificHeatCapacity cp=1;
        parameter SI.Temperature T_ini=293.15;
        parameter Boolean Steadystate_ini=true;
        parameter SI.Temperature T_scale=1;
        parameter Real zeta_a_ini=0;
        parameter Real zeta_b_ini=0;

        Real T(start=T_ini/T_scale) "lumped temperature of the control volume";
        SI.HeatFlowRate Qdot_in;
        SI.HeatFlowRate Qdot_out;
        Real zeta_a(max=1,min=0,start=zeta_a_ini);
        Real zeta_b(max=1,min=0,start=zeta_b_ini);

      protected
        SI.Length L_ab;
        SI.Temperature T_b "Boundary temperature on the right";
        SI.Temperature T_a "Boundary temperature on the left";

      equation

        L_ab=(zeta_b-zeta_a)*L;
      /* Energy balance */
        der(T*T_scale) = (Qdot_in - Qdot_out)/(M/L*cp)/L_ab + (T*T_scale-T_a)/(zeta_b*L)*der(zeta_a*L) +
        (T_b-T*T_scale)/(1-zeta_a)/L*der(zeta_b*L);

      /* Boundary conditions */
        Qdot_in = heatport_in.Q_flow;
        Qdot_out = - heatport_out.Q_flow;
        T*T_scale = heatport_in.T;
        T*T_scale = heatport_out.T;
        heatport_in.L_a=zeta_a*L;
        heatport_out.L_a=zeta_a*L;
        heatport_in.L_b=zeta_b*L;
        heatport_out.L_b=zeta_b*L;

      initial equation
        if Steadystate_ini then
          der(T)=0;
        end if;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end MetalWall;

      partial model IncompressibleLiquid
        replaceable package Medium =
            Modelica.Media.Interfaces.PartialSimpleMedium;                       //Incompressible fluid
        Medium.ThermodynamicState state;
        Medium.ThermodynamicState Inletstate;

        parameter SI.Area HTA=1 "heat transfer are";
        parameter SI.Area Ac=1 "Cross-sectional area";
        parameter SI.Length L=1;
        parameter Medium.Temperature T_ini=Medium.reference_T;
        parameter SI.CoefficientOfHeatTransfer HTC_nominal=1;
        parameter Boolean Steadystate_ini=true;
        parameter Medium.Temperature T_scale=1;
        parameter Real zeta_a_ini=0;
        parameter Real zeta_b_ini=0;

        Medium.AbsolutePressure p;
        parameter Medium.SpecificHeatCapacity cp=Medium.cp_const;
        parameter Medium.Density rho=Medium.d_const;
        Real T(start=T_ini/T_scale);
        Medium.Temperature Tin;
        Medium.SpecificEnthalpy h;
        Medium.Temperature Tt(displayUnit="degC") "Tube temperature";
        Medium.MassFlowRate m_dot;
        Medium.Temperature T_a(displayUnit="degC");
        Medium.Temperature T_b(displayUnit="degC");
        SI.EnergyFlowRate Q_dot;

        SI.Length L_ab;
        Real zeta_a(min=0,max=1,start=zeta_a_ini);
        Real zeta_b(min=0,max=1,start=zeta_b_ini);

      /************************* Connectors *************************/
        TransientVCC.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium);
        TransientVCC.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium);
        DynamicVCC.Components.Units.HX.MB.Interfaces.MB_HeatPort_in heatport;
      equation
        state=Medium.setState_pT(p,T*T_scale);
        h=Medium.specificEnthalpy(state);
        Inletstate=Medium.setState_phX(p,inlet.h_outflow);
        Tin=Medium.temperature(Inletstate);
        //cp=Medium.specificHeatCapacityCp(state);
        //rho=Medium.density(state);

        L_ab=(zeta_b-zeta_a)*L;
        Q_dot=HTA*L_ab/L*HTC_nominal*(Tt-T*T_scale);

      /**************** Mass balance *******************/
        inlet.m_flow+outlet.m_flow=0;

      /**************** Energy balance *******************/
        L_ab*der(T*T_scale)+(T*T_scale-T_b)*der(zeta_b*L)-
        (T*T_scale-T_a)*der(zeta_a*L)=(m_dot*cp*(T_b-T*T_scale)+Q_dot)/(rho*Ac*cp);

      /**************** Boundary Conditions *******************/
        inlet.h_outflow=inStream(inlet.h_outflow);
        m_dot=inlet.m_flow;
        Tt=heatport.T;
        outlet.h_outflow=h;
        heatport.Q_flow=Q_dot;
        inlet.p=p;
        outlet.p=p;
        heatport.L_a=zeta_a*L;
        heatport.L_b=zeta_b*L;

      initial equation
        if Steadystate_ini then
          der(T)=0;
        end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end IncompressibleLiquid;
    end BaseClasses;

    package WorkingFluid
      model SinglePhase
        extends BaseClasses.RefControlVolume;

        /* Variables */
        Medium.Density rho "Average density";
        Medium.Density rho_in;
        Medium.Density rho_out;
        Medium.DerEnthalpyByPressure dhdp_sat "derivative of enthalpy wrt pressure at the saturation line";
        Real dhdt_in "Time derivative of inlet enthalpy";
        Real dhdt_out "Time derivative of outlet enthalpy";
        Real drhodt "Time derivative of average density";

      equation
        rho=Medium.density(state);
        drhodt=drhodp*der(p*p_scale)+drhodh/2*(dhdt_in+dhdt_out);

        if PhaseChange then
          if InletZone then
            if subcooled then
              rho_out=rho_f;
              dhdp_sat=dhdp_f;
              h_out*h_scale=h_f;
            else
              rho_out=rho_g;
              dhdp_sat=dhdp_g;
              h_out*h_scale=h_g;
            end if;
          rho_in=rho;
          dhdt_in=der(h_in);
          dhdt_out=dhdp_sat*der(p*p_scale);
          else
            if subcooled then
              rho_in=rho_f;
              dhdp_sat=dhdp_f;
            else
              rho_in=rho_g;
              dhdp_sat=dhdp_g;
            end if;
            rho_out=rho;
            dhdt_in=dhdp_sat*der(p*p_scale);
            dhdt_out=der(h_out*h_scale);
          end if;
        else //one-zone
          rho_in=rho;
          rho_out=rho;
          dhdp_sat=0;
          dhdt_in=der(h_in);
          dhdt_out=der(h_out*h_scale);
          end if;

      /*
     if InletZone then
      if subcooled then
        rho_out=rho_f;
        dhdp_sat=dhdp_f;
        h_out=h_f;
      else
        rho_out=rho_g;
        dhdp_sat=dhdp_g;
        h_out=h_g;
      end if;
      rho_in=rho;
      dhdt_in=der(h_in);
      dhdt_out=dhdp_sat*der(p*p_scale);
    else
      if subcooled then
        rho_in=rho_f;
        dhdp_sat=dhdp_f;
      else
        rho_in=rho_g;
        dhdp_sat=dhdp_g;
      end if;
      rho_out=rho;
      dhdt_in=dhdp_sat*der(p*p_scale);
      dhdt_out=der(h_out);
    end if;
*/

        /* Mass balance */
        L_ab*drhodt+(rho-rho_out)*der(zeta_b*L)+(rho_in-rho)*der(zeta_a*L)=(m_dot_in-m_dot_out)/Ac;

        /* Energy balance */
        L_ab*(rho*(dhdt_in+dhdt_out)/2+h*drhodt) + (rho*h-rho_out*h_out*h_scale)*der(zeta_b*L) +
        (rho_in*h_in-rho*h)*der(zeta_a*L) - L_ab*der(p*p_scale) =
        (m_dot_in*h_in - m_dot_out*h_out*h_scale - Q_dot)/Ac;

      end SinglePhase;

      model TwoPhase_pVF "Pressure and mean void fraction as states"

        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.RefControlVolume;

        /* Variables */
        Real gamma(min=0,max=1) "Mean void fraction";
        Medium.Density rho_in;
        Medium.Density rho_out;
        parameter Boolean ConstVF=true "Constant mean void fraction";
        parameter Real const_VF=0.8;

      protected
        Medium.DerDensityByPressure drhodp_g;
        Medium.DerDensityByPressure drhodp_f;
        Medium.Density rho_gf;
        Real dgammadt "Time derivative of mean void fraction";
        Real rho_ratio;
        Real x_a "Inlet vapor quality";
        Real x_b "Outlet vapor quality";

      equation

        /* Thermodynamic Properties */
        drhodp_g=Medium.dDewDensity_dPressure(satstate);
        drhodp_f=Medium.dBubbleDensity_dPressure(satstate);
        //rho=gamma*rho_g+(1-gamma)*rho_f;

        rho_gf=rho_g-rho_f;
        rho_ratio=(rho_g/rho_f)^(2/3);
        x_a=(h_in-h_f)/(h_g-h_f);
        x_b=(h_out*h_scale-h_f)/(h_g-h_f);
        dgammadt=0;

        if ConstVF then
          gamma=const_VF;
        else
          gamma=rho_ratio/(x_b-x_a)/(rho_ratio-1)^2*log((x_a*(rho_ratio-1)-rho_ratio)/(x_b*(rho_ratio-1)-rho_ratio))
          -1/(rho_ratio-1);
        end if;

      /*
  if PhaseChange then
    if InletZone then
      if subcooled then
        rho_out=rho_f;
        h_out=h_f;
      else
        rho_out=rho_g;
        h_out=h_g;
      end if;
      rho_in=rho;
    else
      if subcooled then
        rho_in=rho_f;
      else
        rho_in=rho_g;
      end if;
      rho_out=rho;
    end if;
  else
    if InletZone then
      if subcooled then //three-zone condenser
       rho_in=rho_g;
       rho_out=rho_f;
       h_out=h_f;
      else  //three-zone evaporator
        rho_in=rho_f;
        rho_out=rho_g;
        h_out=h_g;
      end if;
    else
      rho_in=rho;
      rho_out=rho;
    end if;
  end if;
  */

          if PhaseChange then
            if InletZone then
              if subcooled then
                rho_out=rho_f;
                h_out*h_scale=h_f;
                rho_in=rho_g;
              else
                rho_out=rho_g;
                h_out*h_scale=h_g;
                rho_in=rho_f;
              end if;
            else
              if subcooled then
                rho_in=rho_f;
              else
                rho_in=rho_g;
              end if;
              rho_out=rho_f;//not used
            end if;
          else
            rho_in=rho_f; //not used
            rho_out=rho_f; //not used
          end if;

        //Mass balance
        L_ab*(gamma*drhodp_g+(1-gamma)*drhodp_f)*der(p*p_scale) +
        L_ab*rho_gf*dgammadt + (rho_f-rho_out+gamma*rho_gf)*der(zeta_b*L) +
        (rho_in-rho_f-rho_gf*gamma)*der(zeta_a*L) = (m_dot_in-m_dot_out)/Ac;

        //Energy balance
        L_ab*der(p*p_scale)*(gamma*(rho_g*dhdp_g+h_g*drhodp_g)+
        (1-gamma)*(rho_f*dhdp_f+h_f*drhodp_f)-1) + L_ab*(rho_g*h_g-rho_f*h_f)*dgammadt +
        (rho_f*h_f-rho_out*h_out+gamma*(rho_g*h_g-rho_f*h_f))*der(L*zeta_b) +
        (rho_in*h_in-rho_f*h_f+gamma*(rho_f*h_f-rho_g*h_g))*der(L*zeta_a) =
        (m_dot_in*h_in-m_dot_out*h_out*h_scale-Q_dot)/Ac;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end TwoPhase_pVF;

      model Inactive

        parameter SI.Length L=1;
        parameter Real zeta_ini=0;
        Real zeta(min=0,max=1,start=zeta_ini,fixed=true);
        parameter Boolean InletZone=true;
        /******************* Connector *******************/
        DynamicVCC.Components.Units.HX.MB.Interfaces.MB_port mb_port;

      equation
        der(zeta)=0;
        zeta*L=mb_port.L;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Inactive;

      model TwoPhase_ph

        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.RefControlVolume;

        /* Variables */
        Real gamma(min=0,max=1) "Mean void fraction";
        Medium.Density rho_in;
        Medium.Density rho_out;
        parameter Boolean ConstVF=true "Constant mean void fraction";
        parameter Real const_VF=0.8;
        Medium.Density rho "Average density";
        Real dhdt_in "Time derivative of inlet enthalpy";
        Real dhdt_out "Time derivative of outlet enthalpy";
        parameter Boolean normal=true; //true if three-zone

      protected
        Medium.DerDensityByPressure drhodp_g;
        Medium.DerDensityByPressure drhodp_f;
        Medium.Density rho_gf;
        Real rho_ratio;
        Real x_a "Inlet vapor quality";
        Real x_b "Outlet vapor quality";
        Real drhodt "Time derivative of average density";
        Real dhdt "Time derivative of average enthalpy";

      equation

        /* Thermodynamic Properties */
        drhodp_g=Medium.dDewDensity_dPressure(satstate);
        drhodp_f=Medium.dBubbleDensity_dPressure(satstate);
        rho=gamma*rho_g+(1-gamma)*rho_f;

        rho_gf=rho_g-rho_f;
        rho_ratio=(rho_g/rho_f)^(2/3);
        x_a=(h_in-h_f)/(h_g-h_f);
        x_b=(h_out*h_scale-h_f)/(h_g-h_f);

        if ConstVF then
          gamma=const_VF;
        else
          gamma=rho_ratio/(x_b-x_a)/(rho_ratio-1)^2*log((x_a*(rho_ratio-1)-rho_ratio)/(x_b*(rho_ratio-1)-rho_ratio))
          -1/(rho_ratio-1);
        end if;

          if PhaseChange then
            if normal then
              if subcooled then //condenser
                rho_out=rho_f;
                h_out*h_scale=h_f;
                rho_in=rho_g;
                dhdt_in=dhdp_g*der(p*p_scale);
                dhdt_out=dhdp_f*der(p*p_scale);
              else             //evaporator
                rho_out=rho_g;
                h_out*h_scale=h_g;
                rho_in=rho_f;
                dhdt_in=dhdp_f*der(p*p_scale);
                dhdt_out=dhdp_g*der(p*p_scale);
              end if;
            else
              if InletZone then
                if subcooled then //TP-SC
                  rho_out=rho_f;
                  h_out*h_scale=h_f;
                  dhdt_out=dhdp_f*der(p*p_scale);
                else
                  rho_out=rho_g;  //TP-SH
                  h_out*h_scale=h_g;
                  dhdt_out=dhdp_g*der(p*p_scale);
                end if;
                dhdt_in=der(h_in);
                rho_in=rho;
              else
                if subcooled then //SC-TP
                  rho_in=rho_f;
                  dhdt_in=dhdp_f*der(p*p_scale);
                else             //SH-TP
                  rho_in=rho_g;
                  dhdt_in=dhdp_g*der(p*p_scale);
                end if;
                rho_out=rho;
                dhdt_out=der(h_out*h_scale);
              end if;
            end if;
          else  //one-zone
            rho_in=rho;
            rho_out=rho;
            dhdt_in=der(h_in);
            dhdt_out=der(h_out*h_scale);
          end if;

        /******************** Governing equations ************************/
        drhodt=drhodp*der(p*p_scale)+drhodh*dhdt;
        dhdt=(dhdt_in+dhdt_out)/2;

        //Mass balance
        L_ab*drhodt+(rho-rho_out)*der(zeta_b*L)+(rho_in-rho)*der(zeta_a*L)=(m_dot_in-m_dot_out)/Ac;

        //Energy balance
        L_ab*h*drhodt+L_ab*rho*dhdt-L_ab*der(p*p_scale) + (rho*h-rho_out*h_out)*der(zeta_b*L) +
         (rho_in*h_in-rho*h)*der(zeta_a*L)=
        (m_dot_in*h_in-m_dot_out*h_out*h_scale-Q_dot)/Ac;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end TwoPhase_ph;
    end WorkingFluid;

    package MetalWall
      model TP
        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.MetalWall;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port_a;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port_b;

      equation
        /* Boundary conditions */
        //zeta_a*L=mb_port_a.L;
        T_a=port_a.T_out;
        //zeta_b*L=mb_port_b.L;
        T_b=port_b.T_out;
        T*T_scale=port_a.T_in;
        T*T_scale=port_b.T_in;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end TP;

      model SP "Tube volume for single phase working fluid volume"
        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.MetalWall;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port;
        parameter Boolean InletZone=true;

      equation

          if InletZone then
            T_a=T*T_scale;
            //zeta_b*L=mb_port.L;
            T_b=port.T_in;
            port.T_out=T*T_scale;
          else
            T_b=T*T_scale;
            //zeta_a*L=mb_port.L;
            T_a=port.T_in;
            port.T_out=T*T_scale;
          end if;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end SP;

      model InactiveZone

        parameter SI.Temperature T_ini=298;
        parameter SI.Temperature T_scale=1;
        Real T(start=T_ini/T_scale);
        SI.Temperature T_b;
        parameter SI.Time tau=1;

        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port;

      equation

            der(T*T_scale)=(T_b-T*T_scale)/tau;
            T*T_scale=port.T_out;
            T_b=port.T_in;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end InactiveZone;
    end MetalWall;

    package IncompressibleLiquid

      model Inlet

        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.IncompressibleLiquid;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port;

      equation

        T_a=port.T_in;
        T_b=Tin;
        port.T_out=T*T_scale;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Inlet;

      model Outlet
        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.IncompressibleLiquid;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port;

      equation

        T_a=T*T_scale;
        T_b=port.T_in;
        port.T_out=T*T_scale;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Outlet;

      model TP
        extends DynamicVCC.Components.Units.HX.MB.Components.BaseClasses.IncompressibleLiquid;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port_a;
        DynamicVCC.Components.Units.HX.MB.Interfaces.T_port port_b;

      equation

        T_b=port_b.T_out;
        T_a=port_a.T_out;
        T*T_scale=port_a.T_in;
        T*T_scale=port_b.T_in;

        //zeta_a*L=mb_port_a.L;
        //zeta_b*L=mb_port_b.L;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end TP;
    end IncompressibleLiquid;
  end Components;

  package Interfaces
    extends Modelica.Icons.InterfacesPackage;
    connector MB_port
      SI.Length L;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end MB_port;

    connector MB_HeatPort_in
      extends TransientVCC.Interfaces.HeatPort_in;
      SI.Length L_a;
      SI.Length L_b;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end MB_HeatPort_in;

    connector MB_HeatPort_out
      extends TransientVCC.Interfaces.HeatPort_out;
      SI.Length L_a;
      SI.Length L_b;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end MB_HeatPort_out;

    connector T_port
      SI.Temperature T_in;
      SI.Temperature T_out;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end T_port;

    model Boundary

      DynamicVCC.Components.Units.HX.MB.Interfaces.MB_port mb_port;
      parameter Boolean InletZone=true;
      SI.Length L=1;

    equation
      if InletZone then
        mb_port.L=0;
      else
        mb_port.L=L;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Boundary;
  end Interfaces;

  package Units
    package Condenser
      partial model SHTPSC "Three-zone representation"
         replaceable package Medium_1 = ExternalMedia.Media.CoolPropMedium constrainedby Modelica.Media.Interfaces.PartialMedium;

        extends DynamicVCC.Components.Units.HX.MB.Units.BaseClass;

         /**************** Components *****************/
        DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.SinglePhase SH(
          redeclare package Medium = Medium_1,
          L=L,
          HTA=HTA,
          Ac=Ac,
          HTC_nominal=HTC_nominal_SP,
          Const_HTC=Const_HTC,
          m_nominal=m_nominal,
          Steadystate_ini=Steadystate_ini,
          p_scale=p_scale,
          p_ini=p_ini,
          zeta_a_ini=0,
          zeta_b_ini=L_a_ini/L,
          PhaseChange=true,
          InletZone=true,
          subcooled=false);

        DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.TwoPhase_ph TP(
          redeclare package Medium = Medium_1,
          L=L,
          HTA=HTA,
          Ac=Ac,
          HTC_nominal=HTC_nominal_TP,
          Const_HTC=Const_HTC,
          m_nominal=m_nominal,
          ConstVF=ConstVF,
          const_VF=const_VF,
          Steadystate_ini=Steadystate_ini,
          p_scale=p_scale,
          p_ini=p_ini,
          PhaseChange=true,
          normal=true,
          subcooled=true,
          zeta_a_ini=L_a_ini/L,
          zeta_b_ini=L_b_ini/L);

        DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.SinglePhase SC(
          redeclare package Medium = Medium_1,
          L=L,
          HTA=HTA,
          Ac=Ac,
          HTC_nominal=HTC_nominal_SP,
          Const_HTC=Const_HTC,
          m_nominal=m_nominal,
          Steadystate_ini=Steadystate_ini,
          p_scale=p_scale,
          h_scale=h_scale,
          p_ini=p_ini,
          h_ini=h_ini,
          zeta_a_ini=L_b_ini/L,
          zeta_b_ini=1,
          PhaseChange=true,
          InletZone=false,
          subcooled=true);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.SP tube_SH(
          M=Mt,
          L=L,
          cp=cp_t,
          T_ini=T_ini_t[1],
          zeta_a_ini=0,
          zeta_b_ini=L_a_ini/L,
          T_scale=T_scale,
          Steadystate_ini=Steadystate_ini,
          InletZone=true);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.TP tube_TP(
          M=Mt,
          L=L,
          cp=cp_t,
          T_ini=T_ini_t[2],
          zeta_a_ini=L_a_ini/L,
          zeta_b_ini=L_b_ini/L,
          T_scale=T_scale,
          Steadystate_ini=Steadystate_ini);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.SP tube_SC(
          M=Mt,
          L=L,
          cp=cp_t,
          T_ini=T_ini_t[3],
          zeta_a_ini=L_b_ini/L,
          zeta_b_ini=1,
          T_scale=T_scale,
          Steadystate_ini=Steadystate_ini,
          InletZone=false);

        DynamicVCC.Components.Units.HX.MB.Components.IncompressibleLiquid.Inlet secondaryfluid_inlet(
          redeclare package Medium = Medium_2,
          HTA=HTA_e,
          Ac=Ac_e,
          L=L,
          T_ini=T_ini_e[1],
          HTC_nominal=HTC_nominal_e,
          Steadystate_ini=Steadystate_ini,
          zeta_a_ini=L_b_ini/L,
          zeta_b_ini=1,
          T_scale=T_scale);

        DynamicVCC.Components.Units.HX.MB.Components.IncompressibleLiquid.TP secondaryfluid_TP(
          redeclare package Medium = Medium_2,
          HTA=HTA_e,
          Ac=Ac_e,
          L=L,
          T_ini=T_ini_e[1],
          HTC_nominal=HTC_nominal_e,
          Steadystate_ini=Steadystate_ini,
          zeta_a_ini=L_a_ini/L,
          zeta_b_ini=L_b_ini/L,
          T_scale=T_scale);

        DynamicVCC.Components.Units.HX.MB.Components.IncompressibleLiquid.Outlet secondaryfluid_outlet(
          redeclare package Medium = Medium_2,
          HTA=HTA_e,
          Ac=Ac_e,
          L=L,
          T_ini=T_ini_e[1],
          HTC_nominal=HTC_nominal_e,
          Steadystate_ini=Steadystate_ini,
          zeta_a_ini=0,
          zeta_b_ini=L_a_ini/L,
          T_scale=T_scale);

        SI.EnergyFlowRate Qtot_dot;

      equation

        connect(secondaryfluid_outlet.heatport,tube_SH.heatport_out);
        connect(secondaryfluid_TP.heatport,tube_TP.heatport_out);
        connect(secondaryfluid_inlet.heatport,tube_SC.heatport_out);
        connect(secondaryfluid_outlet.port,secondaryfluid_TP.port_a);
        connect(secondaryfluid_TP.port_b,secondaryfluid_inlet.port);

        connect(secondaryfluid_inlet.outlet,secondaryfluid_TP.inlet);
        connect(secondaryfluid_TP.outlet,secondaryfluid_outlet.inlet);

        connect(SH.outlet,TP.inlet);
        connect(TP.outlet,SC.inlet);
        connect(boundary_a.mb_port,SH.mb_port_a);
        connect(SH.mb_port_b,TP.mb_port_a);
        connect(TP.mb_port_b,SC.mb_port_a);
        connect(SC.mb_port_b,boundary_b.mb_port);

        connect(SH.heatport,tube_SH.heatport_in);
        connect(TP.heatport,tube_TP.heatport_in);
        connect(SC.heatport,tube_SC.heatport_in);
        connect(tube_SH.port,tube_TP.port_a);
        connect(tube_TP.port_b,tube_SC.port);

        Qtot_dot=abs(SH.Q_dot+TP.Q_dot+SC.Q_dot);

      initial equation
        if Steadystate_ini then
          der(SC.h_out)=0;
        end if;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end SHTPSC;
    end Condenser;

    package Evaporator
      model TPSH

        extends DynamicVCC.Components.Units.HX.MB.Units.BaseClass;
         /**************** Components *****************/

        DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.TwoPhase_ph TP(
          redeclare package Medium = Medium_1,
          L=L,
          HTA=HTA,
          Ac=Ac,
          HTC_nominal=HTC_nominal_TP,
          Const_HTC=Const_HTC,
          m_nominal=m_nominal,
          ConstVF=ConstVF,
          const_VF=const_VF,
          Steadystate_ini=Steadystate_ini,
          p_scale=p_scale,
          p_ini=p_ini,
          PhaseChange=true,
          normal=false,
          InletZone=true,
          subcooled=false,
          zeta_a_ini=0,
          zeta_b_ini=L_b_ini/L);

        DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.SinglePhase SH(
          redeclare package Medium = Medium_1,
          L=L,
          HTA=HTA,
          Ac=Ac,
          HTC_nominal=HTC_nominal_SP,
          Const_HTC=Const_HTC,
          m_nominal=m_nominal,
          Steadystate_ini=Steadystate_ini,
          p_scale=p_scale,
          p_ini=p_ini,
          h_ini=h_ini,
          zeta_a_ini=L_b_ini/L,
          zeta_b_ini=1,
          PhaseChange=true,
          InletZone=false,
          subcooled=false);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.InactiveZone tube_SC(T_scale=T_scale, T_ini=T_ini_t[1]);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.TP tube_TP(
          M=Mt,
          L=L,
          cp=cp_t,
          T_ini=T_ini_t[2],
          zeta_a_ini=0,
          zeta_b_ini=L_b_ini/L,
          T_scale=T_scale,
          Steadystate_ini=Steadystate_ini);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.SP tube_SH(
          M=Mt,
          L=L,
          cp=cp_t,
          T_ini=T_ini_t[3],
          zeta_a_ini=L_b_ini/L,
          zeta_b_ini=1,
          T_scale=T_scale,
          Steadystate_ini=Steadystate_ini,
          InletZone=false);

        DynamicVCC.Components.Units.HX.MB.Components.IncompressibleLiquid.Inlet secondaryfluid_inlet(
          redeclare package Medium = Medium_2,
          HTA=HTA_e,
          Ac=Ac_e,
          L=L,
          T_ini=T_ini_e[1],
          HTC_nominal=HTC_nominal_e,
          Steadystate_ini=Steadystate_ini,
          T_scale=T_scale);

        DynamicVCC.Components.Units.HX.MB.Components.IncompressibleLiquid.TP secondaryfluid_TP(
          redeclare package Medium = Medium_2,
          HTA=HTA_e,
          Ac=Ac_e,
          L=L,
          T_ini=T_ini_e[2],
          HTC_nominal=HTC_nominal_e,
          Steadystate_ini=Steadystate_ini,
          T_scale=T_scale);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.InactiveZone secondaryfluid_outlet(T_scale=T_scale, T_ini=T_ini_e[3]);

        TransientVCC.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium_1);
        TransientVCC.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium_1);
        TransientVCC.Interfaces.FluidPort_a secondary_inlet(redeclare package Medium = Medium_2);
        TransientVCC.Interfaces.FluidPort_b secondary_outlet(redeclare package Medium = Medium_2);

      equation
        connect(TP.outlet,SH.inlet);
        connect(boundary_a.mb_port,TP.mb_port_a);
        connect(TP.mb_port_b,SH.mb_port_a);
        connect(SH.mb_port_b,boundary_b.mb_port);

        connect(TP.heatport,tube_TP.heatport_in);
        connect(SH.heatport,tube_SH.heatport_in);
        connect(tube_SC.port,tube_TP.port_a);
        connect(tube_TP.port_b,tube_SH.port);

        connect(secondaryfluid_TP.heatport,tube_TP.heatport_out);
        connect(secondaryfluid_inlet.heatport,tube_SH.heatport_out);
        connect(secondaryfluid_outlet.port,secondaryfluid_TP.port_a);
        connect(secondaryfluid_TP.port_b,secondaryfluid_inlet.port);
        connect(secondaryfluid_inlet.outlet,secondaryfluid_TP.inlet);

        connect(secondary_inlet,secondaryfluid_inlet.inlet);
        connect(secondaryfluid_TP.outlet,secondary_outlet);
        connect(inlet,TP.inlet);
        connect(SH.outlet,outlet);

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end TPSH;

      model TP

        extends DynamicVCC.Components.Units.HX.MB.Units.BaseClass;
         /**************** Components *****************/

        DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.TwoPhase_ph TP(
          redeclare package Medium = Medium_1,
          L=L,
          HTA=HTA,
          Ac=Ac,
          HTC_nominal=HTC_nominal_TP,
          Const_HTC=Const_HTC,
          m_nominal=m_nominal,
          ConstVF=ConstVF,
          const_VF=const_VF,
          Steadystate_ini=Steadystate_ini,
          p_scale=p_scale,
          h_scale=h_scale,
          p_ini=p_ini,
          PhaseChange=false,
          zeta_a_ini=0,
          zeta_b_ini=1);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.InactiveZone tube_SC(T_scale=T_scale, T_ini=T_ini_t[1]);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.TP tube_TP(
          M=Mt,
          L=L,
          cp=cp_t,
          T_ini=T_ini_t[2],
          zeta_a_ini=0,
          zeta_b_ini=L_b_ini/L,
          T_scale=T_scale,
          Steadystate_ini=Steadystate_ini);

      DynamicVCC.Components.Units.HX.MB.Components.MetalWall.InactiveZone tube_SH(T_scale=T_scale, T_ini=T_ini_t[3]);

          DynamicVCC.Components.Units.HX.MB.Components.MetalWall.InactiveZone secondaryfluid_inlet(T_scale=T_scale, T_ini=T_ini_e[1]);

        DynamicVCC.Components.Units.HX.MB.Components.IncompressibleLiquid.TP secondaryfluid_TP(
          redeclare package Medium = Medium_2,
          HTA=HTA_e,
          Ac=Ac_e,
          L=L,
          T_ini=T_ini_e[2],
          HTC_nominal=HTC_nominal_e,
          Steadystate_ini=Steadystate_ini,
          T_scale=T_scale);

        DynamicVCC.Components.Units.HX.MB.Components.MetalWall.InactiveZone secondaryfluid_outlet(T_scale=T_scale, T_ini=T_ini_e[3]);

        TransientVCC.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium_1);
        TransientVCC.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium_1);
        TransientVCC.Interfaces.FluidPort_a secondary_inlet(redeclare package Medium = Medium_2);
        TransientVCC.Interfaces.FluidPort_b secondary_outlet(redeclare package Medium = Medium_2);

      equation

        connect(boundary_a.mb_port,TP.mb_port_a);
        connect(TP.mb_port_b,boundary_b.mb_port);
        connect(inlet,TP.inlet);
        connect(TP.outlet,outlet);

        connect(TP.heatport,tube_TP.heatport_in);
        connect(tube_SC.port,tube_TP.port_a);
        connect(tube_TP.port_b,tube_SH.port);

        connect(secondaryfluid_TP.heatport,tube_TP.heatport_out);

        connect(secondary_inlet,secondaryfluid_inlet.inlet);
        connect(secondaryfluid_TP.inlet,secondaryfluid_inlet.outlet);
        connect(secondaryfluid_TP.outlet,secondaryfluid_outlet.inlet);
        connect(secondaryfluid_outlet.outlet,secondary_outlet);

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end TP;
    end Evaporator;

    partial model BaseClass
        replaceable package Medium_1 = ExternalMedia.Media.CoolPropMedium constrainedby Modelica.Media.Interfaces.PartialMedium;

       /**************** Geometry ****************/
       parameter SI.Length L=1;
       parameter SI.Area HTA=1 "Total heat transfer are of working fluid";
       parameter SI.Area Ac=1;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_TP=1;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_SP=1;
       parameter Boolean Steadystate_ini=true;
       parameter Boolean ConstVF=true "Constant mean void fraction";
       parameter Real const_VF=0.8;
       parameter Boolean Const_HTC=true;
       parameter Medium_1.MassFlowRate m_nominal=1;

       parameter SI.Mass Mt=1 "Total mass of tube";
       parameter SI.SpecificHeatCapacity cp_t=1 "Tube material specific heat";

       /*************** Initial conditions *******************/
       parameter Medium_1.AbsolutePressure p_ini=Medium_1.reference_p;
       parameter SI.Temperature T_ini_t[3]=Medium_1.reference_T*ones(3);
       parameter Medium_1.SpecificEnthalpy h_ini=Medium_1.h_default;
       parameter Medium_1.AbsolutePressure p_scale=1;
       parameter Medium_1.SpecificEnthalpy h_scale=1;
       parameter SI.Temperature T_scale=1;
       parameter SI.Length L_a_ini;
       parameter SI.Length L_b_ini;

       /***************** Secondary fluid **********************/
      replaceable package Medium_2 =
          Modelica.Media.Interfaces.PartialSimpleMedium;
      parameter SI.Area HTA_e=1;
      parameter SI.Area Ac_e=1;
      parameter Medium_2.Temperature T_ini_e[3]=Medium_2.reference_T*ones(3);
      parameter SI.CoefficientOfHeatTransfer HTC_nominal_e=1;

      DynamicVCC.Components.Units.HX.MB.Interfaces.Boundary boundary_a(InletZone=true, L=L);
      DynamicVCC.Components.Units.HX.MB.Interfaces.Boundary boundary_b(InletZone=false, L=L);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end BaseClass;
  end Units;

  package Tests
    extends Modelica.Icons.ExamplesPackage;

    model Condenser "Counter flow"

      extends DynamicVCC.Components.Units.HX.MB.Units.Condenser.SHTPSC;

      /******************* Connector ***********************/
      TransientVCC.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium_1);
      TransientVCC.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium_1);
      TransientVCC.Interfaces.FluidPort_a secondary_inlet(redeclare package Medium = Medium_2);
      TransientVCC.Interfaces.FluidPort_b secondary_outlet(redeclare package Medium = Medium_2);

    equation

      connect(inlet,SH.inlet);
      connect(SC.outlet,outlet);
      connect(secondary_inlet,secondaryfluid_inlet.inlet);
      connect(secondaryfluid_outlet.outlet,secondary_outlet);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Condenser;

    model testhx

    replaceable package Medium_1 = ExternalMedia.Media.CoolPropMedium constrainedby Modelica.Media.Interfaces.PartialMedium;

       /**************** Geometry ****************/
       parameter SI.Length L=1;
       parameter SI.Area HTA=1 "Total heat transfer are of working fluid";
       parameter SI.Area Ac=1;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal=1;
       parameter Boolean Steadystate_ini=true;
       parameter Integer N=3 "Active zones";

       parameter SI.Mass Mt=1 "Total mass of tube";
       parameter SI.SpecificHeatCapacity cp_t=1 "Tube material specific heat";

       /*************** Initial conditions *******************/
       parameter Medium_1.AbsolutePressure p_ini=Medium_1.reference_p;
       parameter SI.Temperature Tt_ini[N]=Medium_1.reference_T*ones(N);
       parameter Medium_1.SpecificEnthalpy h_ini=Medium_1.h_default;
       parameter Medium_1.AbsolutePressure p_scale=1;
       parameter Medium_1.SpecificEnthalpy h_scale=1;
       parameter SI.Temperature T_scale=1;
       parameter SI.Length L_a_ini;
       parameter SI.Length L_b_ini;
      /***************** Secondary fluid **********************/
      replaceable package Medium_2 =
          Modelica.Media.Interfaces.PartialSimpleMedium;
      parameter SI.Area HTA_e=1;
      parameter SI.Area Ac_e=1;
      parameter Medium_2.Temperature T_ini_e[N]=Medium_2.reference_T*ones(N);
      parameter SI.CoefficientOfHeatTransfer HTC_nominal_e=1;

      DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.SinglePhase SH(
        redeclare package Medium = Medium_1,
        L=L,
        HTA=HTA,
        Ac=Ac,
        HTC_nominal=HTC_nominal,
        Steadystate_ini=Steadystate_ini,
        p_scale=p_scale,
        p_ini=p_ini,
        zeta_a_ini=0,
        zeta_b_ini=L_a_ini/L,
        InletZone=true,
        subcooled=false);

      DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.TwoPhase_pVF TP(
        redeclare package Medium = Medium_1,
        L=L,
        HTA=HTA,
        Ac=Ac,
        HTC_nominal=HTC_nominal,
        Steadystate_ini=Steadystate_ini,
        p_scale=p_scale,
        p_ini=p_ini,
        PhaseChange=false,
        InletZone=true,
        zeta_a_ini=L_a_ini/L,
        zeta_b_ini=L_b_ini/L);

       DynamicVCC.Components.Units.HX.MB.Components.WorkingFluid.SinglePhase SC(
        redeclare package Medium = Medium_1,
        L=L,
        HTA=HTA,
        Ac=Ac,
        HTC_nominal=HTC_nominal,
        Steadystate_ini=Steadystate_ini,
        p_scale=p_scale,
        p_ini=p_ini,
        h_ini=h_ini,
        zeta_a_ini=L_b_ini/L,
        zeta_b_ini=1,
        InletZone=false,
        subcooled=true);

    DynamicVCC.Components.Units.HX.MB.Tests.HeatSource_T heatsource_SH;
    DynamicVCC.Components.Units.HX.MB.Tests.HeatSource_T heatsource_TP;
    DynamicVCC.Components.Units.HX.MB.Tests.HeatSource_T heatsource_SC;

    DynamicVCC.Components.Units.HX.MB.Interfaces.Boundary boundary_a(InletZone=true);
    DynamicVCC.Components.Units.HX.MB.Interfaces.Boundary boundary_b(InletZone=false);

      /******************* Connector ***********************/
      TransientVCC.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium_1);
      TransientVCC.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium_1);

      TransientVCC.Component.Sources.MassFlowSource_h ref_source;

    equation

      connect(SH.heatport,heatsource_SH.heatport);
      connect(TP.heatport,heatsource_TP.heatport);
      connect(SC.heatport,heatsource_SC.heatport);

      connect(inlet,SH.inlet);
      connect(SH.outlet,TP.inlet);
      connect(TP.outlet,SC.inlet);
      connect(outlet,SC.outlet);
      connect(boundary_a.mb_port,SH.mb_port_a);
      connect(SH.mb_port_b,TP.mb_port_a);
      connect(TP.mb_port_b,SC.mb_port_a);
      connect(boundary_b.mb_port,SC.mb_port_b);
      connect(ref_source.port,inlet);
      //ref_source.h_in=5e5;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end testhx;

    model HeatSource_T
      DynamicVCC.Components.Units.HX.MB.Interfaces.MB_HeatPort_in heatport;
      Modelica.Blocks.Interfaces.RealInput T_in;
      //SI.Length La;
      //SI.Length Lb;
    equation
      heatport.T=T_in;
      //La=heatport.L_a;
      //Lb=heatport.L_b;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end HeatSource_T;

    model Testcondenser
      extends Modelica.Icons.Example;

      package R134a "R134a from CoolProp"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName="R134a",
          substanceNames={
              "R134a|rho_smoothing_xend=0.1|calc_transport=0|enable_TTSE=1"},
          ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
      end R134a;

      package Medium_1=R134a;
      package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater (
      cp_const=4186.8,
      d_const=995);

       /**************** Geometry ****************/
       parameter SI.Length L=2.4384;
       parameter SI.Area HTA_r=23.9328 "Total heat transfer are of working fluid";
       parameter SI.Area Ac_r=0.261;
       parameter Boolean Steadystate_ini=true;
       parameter Boolean ConstVF=false "Constant mean void fraction";
       parameter Real const_VF=0.98;
       parameter Boolean Const_HTC=false;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_TP=0.8e5;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_SP=0.2e5;
       parameter Medium_1.MassFlowRate m_nominal=1.5;
       parameter SI.Mass Mt=340.638 "Total mass of tube";
       parameter SI.SpecificHeatCapacity cp_t=385 "Tube material specific heat";
       parameter SI.Area HTA_e=22.3716;
       parameter SI.Area Ac_e=0.305;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_e=0.3e5;

       /*************** Initial conditions *******************/
       parameter Medium_1.AbsolutePressure p_ini=9.2386e5;
       parameter SI.Temperature T_ini_t[3]={311.1485,309.3957,303.7665};
       parameter Medium_1.SpecificEnthalpy h_ini=2.3671e5 "Exit enthalpy";
       parameter Medium_1.AbsolutePressure p_scale=Medium_1.reference_p;
       parameter Medium_1.SpecificEnthalpy h_scale=Medium_1.h_default;
       parameter SI.Temperature T_scale=Medium_1.reference_T;
       parameter SI.Length L_a_ini=0.0146*L;
       parameter SI.Length L_b_ini=0.9212*L;
       parameter Medium_2.Temperature T_ini_e[3]={303.4388,309.0674,309.5181};

       DynamicVCC.Components.Units.HX.MB.Tests.Condenser condenser(
        redeclare package Medium_1 = Medium_1,
        redeclare package Medium_2 = Medium_2,
        L=L,
        HTA=HTA_r,
        Ac=Ac_r,
        HTC_nominal_SP=HTC_nominal_SP,
        HTC_nominal_TP=HTC_nominal_TP,
        Const_HTC=Const_HTC,
        m_nominal=m_nominal,
        Steadystate_ini=Steadystate_ini,
        p_ini=p_ini,
        h_ini=h_ini,
        ConstVF=ConstVF,
        const_VF=const_VF,
        T_ini_t=T_ini_t,
        Mt=Mt,
        cp_t=cp_t,
        HTA_e=HTA_e,
        Ac_e=Ac_e,
        T_ini_e=T_ini_e,
        HTC_nominal_e=HTC_nominal_e,
        L_a_ini=L_a_ini,
        L_b_ini=L_b_ini,
        p_scale=p_scale,
        h_scale=h_scale,
        T_scale=T_scale);

      Modelica.Blocks.Sources.CombiTimeTable Cond_BC(tableOnFile=true,smoothness=5,tableName="Cond_BC",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Cond_BC.mat",columns=2:6);
      Modelica.Blocks.Sources.CombiTimeTable Mea_Cond(tableOnFile=true,smoothness=1,tableName="Mea_Cond",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Mea_Cond.mat",columns=2:5);
      Modelica.Blocks.Sources.CombiTimeTable Evap_BC(tableOnFile=true,smoothness=1,tableName="Evap_BC",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Evap_BC.mat",columns=2:6);
      Modelica.Blocks.Sources.CombiTimeTable Mea_Evap(tableOnFile=true,smoothness=1,tableName="Mea_Evap",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Mea_Evap.mat",columns=2:5);

     TransientVCC.Component.Sources.Fluid.MassFlowSource_hX ref_in(redeclare package Medium=Medium_1,
      nPorts=1,
      use_m_flow_in=true,
      use_h_in=true);
      TransientVCC.Component.Sources.Fluid.MassFlowSource_hX ref_out(redeclare package Medium=Medium_1,
      nPorts=1,
      use_m_flow_in=true);
      TransientVCC.Component.Sources.Fluid.MassFlowSource_TX water_in(redeclare package Medium=Medium_2,
      nPorts=1,
      use_m_flow_in=true,
      use_T_in=true);
      Modelica.Fluid.Sources.Boundary_pT water_out(redeclare package Medium = Medium_2,
      nPorts=1);

    equation

      water_in.T_in=Cond_BC.y[4];
      water_in.m_flow_in=16.95;
      ref_in.h_in=Cond_BC.y[3];
      ref_in.m_flow_in=Cond_BC.y[1];
      ref_out.m_flow_in=-Cond_BC.y[2];

    /*
  water_in.T_in=302.95;
  water_in.m_flow_in=16.95;
  ref_in.h_in=431209.6;
  ref_in.m_flow_in=2.3964;
  ref_out.m_flow_in=2.3964;
*/
      connect(ref_in.ports[1],condenser.inlet);
      connect(condenser.outlet,ref_out.ports[1]);
      connect(condenser.secondary_inlet,water_in.ports[1]);
      connect(condenser.secondary_outlet,water_out.ports[1]);

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Testcondenser;

    model Evaporator

      extends DynamicVCC.Components.Units.HX.MB.Units.Evaporator.TPSH;

      /******************* Connector ***********************/

    equation

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Evaporator;

    model Testevaporator
      extends Modelica.Icons.Example;

       package R134a "R134a from CoolProp"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName="R134a",
          substanceNames={
              "R134a|rho_smoothing_xend=0.1|calc_transport=0|enable_TTSE=1"},
          ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
       end R134a;

      package Medium_1=R134a;
      package Medium_2=Modelica.Media.Water.ConstantPropertyLiquidWater (
      cp_const=4186.8,
      d_const=995);

       /**************** Geometry ****************/
       parameter SI.Length L=2.4384;
       parameter SI.Area HTA_r=22.3716 "Total heat transfer are of working fluid";
       parameter SI.Area Ac_r=0.0763;
       parameter Boolean Steadystate_ini=true;
       parameter Boolean ConstVF=false "Constant mean void fraction";
       parameter Real const_VF=0.98;
       parameter Boolean Const_HTC=true;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_TP=0.8e5;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_SP=0.2e5;
       parameter Medium_1.MassFlowRate m_nominal=1.5;
       parameter SI.Mass Mt=321.783 "Total mass of tube";
       parameter SI.SpecificHeatCapacity cp_t=385 "Tube material specific heat";
       parameter SI.Area HTA_e=18.3310;
       parameter SI.Area Ac_e=0.045;
       parameter SI.CoefficientOfHeatTransfer HTC_nominal_e=0.4e5;

       /*************** Initial conditions *******************/
       parameter Medium_1.AbsolutePressure p_ini=3.9e5;
       parameter SI.Temperature T_ini_t[3]={282,284,286};
       parameter Medium_1.SpecificEnthalpy h_ini=406e3 "Exit enthalpy";
       parameter Medium_1.AbsolutePressure p_scale=Medium_1.p_default;
       parameter Medium_1.SpecificEnthalpy h_scale=1;
       parameter SI.Temperature T_scale=Medium_1.T_default;
       parameter SI.Length L_a_ini=0*L;
       parameter SI.Length L_b_ini=0.99*L;
       parameter Medium_2.Temperature T_ini_e[3]={288,285,283};

       //TransientVCC.Component.Units.HX.MB.Tests.Evaporator evaporator(
       DynamicVCC.Components.Units.HX.MB.Units.Evaporator.TPSH evaporator(
        redeclare package Medium_1 = Medium_1,
        redeclare package Medium_2 = Medium_2,
        L=L,
        HTA=HTA_r,
        Ac=Ac_r,
        HTC_nominal_SP=HTC_nominal_SP,
        HTC_nominal_TP=HTC_nominal_TP,
        Const_HTC=Const_HTC,
        m_nominal=m_nominal,
        Steadystate_ini=Steadystate_ini,
        p_ini=p_ini,
        h_ini=h_ini,
        ConstVF=ConstVF,
        const_VF=const_VF,
        T_ini_t=T_ini_t,
        Mt=Mt,
        cp_t=cp_t,
        HTA_e=HTA_e,
        Ac_e=Ac_e,
        T_ini_e=T_ini_e,
        HTC_nominal_e=HTC_nominal_e,
        L_a_ini=L_a_ini,
        L_b_ini=L_b_ini);

      Modelica.Blocks.Sources.CombiTimeTable Cond_BC(tableOnFile=true,smoothness=5,tableName="Cond_BC",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Cond_BC.mat",columns=2:6);
      Modelica.Blocks.Sources.CombiTimeTable Mea_Cond(tableOnFile=true,smoothness=1,tableName="Mea_Cond",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Mea_Cond.mat",columns=2:5);
      Modelica.Blocks.Sources.CombiTimeTable Evap_BC(tableOnFile=true,smoothness=5,tableName="Evap_BC",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Evap_BC.mat",columns=2:6);
      Modelica.Blocks.Sources.CombiTimeTable Mea_Evap(tableOnFile=true,smoothness=1,tableName="Mea_Evap",fileName="C:/Users/alaverni/Documents/Dymola 2019/Jiacheng/Jiacheng_Ma/TPWL/BC/chiller/Mea_Evap.mat",columns=2:5);

      TransientVCC.Component.Sources.Fluid.MassFlowSource_hX ref_in(redeclare package Medium=Medium_1,
      nPorts=1,
      use_m_flow_in=true,
      use_h_in=true);
      TransientVCC.Component.Sources.Fluid.MassFlowSource_hX ref_out(redeclare package Medium=Medium_1,
      nPorts=1,
      use_m_flow_in=true);
      TransientVCC.Component.Sources.Fluid.MassFlowSource_TX water_in(redeclare package Medium=Medium_2,
      nPorts=1,
      use_m_flow_in=true,
      use_T_in=true);
      Modelica.Fluid.Sources.Boundary_pT water_out(redeclare package Medium = Medium_2,
      nPorts=1);

    equation

      water_in.T_in=Evap_BC.y[4];
      water_in.m_flow_in=13.565;
      ref_in.h_in=Evap_BC.y[3];
      ref_in.m_flow_in=Evap_BC.y[1];
      ref_out.m_flow_in=-Evap_BC.y[2];

    /*
  //water_in.T_in=302.95;
  //water_in.m_flow_in=16.95;
  ref_in.h_in=431209.6;
  ref_in.m_flow_in=2.3964;
  ref_out.m_flow_in=2.3964;
*/
      connect(ref_in.ports[1],evaporator.inlet);
      connect(evaporator.outlet,ref_out.ports[1]);
      connect(evaporator.secondary_inlet,water_in.ports[1]);
      connect(evaporator.secondary_outlet,water_out.ports[1]);

    end Testevaporator;
  end Tests;
  annotation (Icon(graphics={
        Rectangle(
          extent={{-80,-42},{-40,20}},
          lineColor={28,108,200},
          fillColor={238,46,47},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-40,20},{40,-42}},
          lineColor={0,255,255},
          fillColor={0,128,255},
          fillPattern=FillPattern.Forward),
        Rectangle(
          extent={{40,-42},{80,20}},
          lineColor={28,108,200},
          fillColor={102,44,145},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-40,20},{-12,6},{18,-4},{50,-8},{6,-20},{-20,-30},{-40,-42}},
          color={238,46,47},
          smooth=Smooth.Bezier,
          thickness=0.5)}));

end MB;
