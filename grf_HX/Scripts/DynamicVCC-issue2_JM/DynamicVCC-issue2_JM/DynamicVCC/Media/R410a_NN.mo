within DynamicVCC.Media;
package R410a_NN "Neural network models for calculating R410A properties"

    extends DynamicVCC.Media.BaseClasses(
    ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph,
    mediumName="R410a",
    substanceNames={"R410a"},
    singleState=false,
    SpecificEnthalpy(
      start=h_default,
      nominal=5e5,
      min=1e5,
      max=5e5),
    Density(start=500, nominal=500),
    AbsolutePressure(
      start=p_default,
      nominal=40e5,
      max=48e5,
      min=1e5),
    Temperature(start=T_default, nominal=350),
    smoothModel=false,
    onePhase=false,
    reference_p=p_default,
    reference_T=T_default,
    p_default=10e5,
    T_default=298,
    h_default=4e5,
    fluidConstants=r410aConstants);

    constant TwoPhase.FluidConstants[1] r410aConstants(
    each iupacName="",
    each casRegistryNumber="",
    each chemicalFormula="",
    each structureFormula="",
    each molarMass=0.0725854,
    each criticalTemperature=344.494,
    each criticalPressure=4901200,
    each criticalMolarVolume=1.1478e-05,
    each acentricFactor=0.296,
    each triplePointTemperature=200,
    each triplePointPressure=29160.3353747,
    each meltingPoint=118.15,
    each normalBoilingPoint=224.65,
    each dipoleMoment=0,
    each hasCriticalData=true);

    extends Modelica.Icons.VariantsPackage;

  redeclare function extends saturationTemperature
    "Mixture has a varying two-phase temperature. For R410A the dew point temperature is returned"
  algorithm
    T:=0.00001713511163982163*p - 5.482231475859701*(0.0000007481855328428306*p - 1.859241049114434)^2 - 0.1085719864892459*(0.0000007481855328428306*p - 1.859241049114434)^3 + 0.804769360499032*(0.0000007481855328428306*p - 1.859241049114434)^4 + 0.9796247613048392*(0.0000007481855328428306*p - 1.859241049114434)^5 - 0.608342334436711*(0.0000007481855328428306*p - 1.859241049114434)^6 + 271.773473897122;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end saturationTemperature;

  redeclare replaceable model extends BaseProperties

    FixedPhase phase;
  equation
    state=setState_ph(p,h);
    sat=setSat_p(p);
    MM=0.07258540000000001;
    R_s=114.55;
    d=state.d;
    T=state.T;
    u=h-p/d;
    phase=state.phase;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end BaseProperties;

  redeclare function extends bubbleEnthalpy
  protected
    parameter AbsolutePressure pc=4901200;
    Real p;
  algorithm
    p:=sat.psat/pc;
    hl:=4e5*(8.8106528/(exp(-33.109335*p - 4.4749095) + 1.0) - 8.2356376/(exp(0.2599934*p - 1.5963615) + 1.0) - 8.7313665/(exp(19.12902*p - 24.506751) + 1.0) - 9.557574/(exp(5.7577211*p + 3.9564407) + 1.0) + 7.2900948)
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));

  end bubbleEnthalpy;

  redeclare function extends dewEnthalpy
  protected
    parameter AbsolutePressure pc=4901200;
    Real x;
  algorithm
    x:=sat.psat/pc;
    hv:=4e5 * (1.1024442/(exp(- 6.9431234*x - 2.6964151) + 1.0) - 0.99955293/(exp(37.48401*x + 3.1660664) + 1.0) - 1.071489/(exp(29.951489 - 26.44918*x) + 1.0) - 1.47877/(exp(6.5300051 - 3.9647703*x) + 1.0) - 0.021353525)
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));

  end dewEnthalpy;

  redeclare function extends bubbleDensity
  protected
    parameter AbsolutePressure pc=4901200;
    Real x;
  algorithm
    x:=sat.psat/pc;
    dl:=459.03*(47.254718 - 23.606547/(exp(- 2.4866747*x - 3.3233854) + 1.0) - 23.134726/(exp(6.1600514 - 2.741962*x) + 1.0) - 21.541987/(exp(30.113038 - 25.591383*x) + 1.0) - 21.580232/(exp(- 24.855288*x - 4.5370323) + 1.0));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end bubbleDensity;

  redeclare function extends dewDensity
  protected
    parameter AbsolutePressure pc=4901200;
    Real x;
  algorithm
    x:=sat.psat/pc;
    dv:=459.03*(45.59568 - 14.32777/(exp(0.62756909*x - 3.109196) + 1.0) - 14.670827/(exp(29.418339*x - 34.105511) + 1.0) - 17.207798/(exp(5.8014555*x - 9.9767952) + 1.0));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end dewDensity;

  redeclare function extends setDewState
  algorithm
    state.p:=sat.psat;
    state.T:=saturationTemperature(sat.psat);
    state.d:=dewDensity(sat);
    state.h:=dewEnthalpy(sat);
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end setDewState;

  redeclare function extends setBubbleState
  protected
    Real p=sat.psat;
  algorithm
    state.p:=sat.psat;
    state.h:=bubbleEnthalpy(sat);
    state.d:=bubbleDensity(sat);
    state.T:=378.3805941655512819928829353219 - 161.45848521884219762203836331733/(exp(0.00000034517364925262776894808201071591*p - 0.36215572886618220142499254810919) + 1.0) - 200.38823078782339219641833959073/(exp(0.0000043736390471608007228016744337754*p + 1.870447838042673240894920055337) + 1.0) - 223.39746953272946690646669568526/(exp(0.0000011561072372861210781043256818345*p + 1.38807816673089497848261793435) + 1.0);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end setBubbleState;

  redeclare function extends specificHeatCapacityCp
   import DynamicVCC.Media.Utilities.phaseTransition;
  protected
    parameter AbsolutePressure pc=4901200;
    parameter SpecificEnthalpy h_ref=4e5;
    AbsolutePressure pr=p/pc;
    SpecificEnthalpy hr=h/h_ref;
  algorithm
    (w1,w2,w3):=phaseTransition(x);
    cp_f:=27411492.0/(exp(138.36401 - 133.3942*pr) + 1.0) - 7301310.4/(exp(12.244959*pr - 18.984874) + 1.0) - 2738399.1/(exp(1.6307769*pr - 8.7149415) + 1.0) + 12199426.0/(exp(50.04042 - 44.300356*pr) + 1.0) + 10040634.0;
    cp_g:=524730.04/(exp(-0.036034187*pr - 2.1569229) + 1.0) - 276914.58/(exp(8.2058641*pr - 11.477291) + 1.0) + 451909.15/(exp(34.510777 - 31.939083*pr) + 1.0) - 759303.31/(exp(104.9862*pr - 106.89429) + 1.0) + 566716.77;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end specificHeatCapacityCp;

  redeclare function extends dynamicViscosity
   import DynamicVCC.Media.Utilities.phaseTransition;
  protected
    parameter AbsolutePressure pc=4901200;
    parameter SpecificEnthalpy h_ref=4e5;
    AbsolutePressure pr=p/pc;
    SpecificEnthalpy hr=h/h_ref;
  algorithm
    (w1,w2,w3):=phaseTransition(x);
    eta_f:=0.044155251/(exp(11.639973*pr - 19.631438) + 1.0) - 0.030698551/(exp(-10.511648*pr - 5.6233328) + 1.0) - 0.020072014/(exp(-1.6604206*pr - 4.8840398) + 1.0) + 0.03260374/(exp(46.294695*pr + 5.7591829) + 1.0) + 0.0066450215;
    eta_g:=0.0041205319/(exp(11.24613 - 5.2345529*pr) + 1.0) - 0.0028436496/(exp(1.1195385*pr + 5.6685719) + 1.0) - 0.003611101/(exp(29.676619*pr - 36.343012) + 1.0) + 0.0034287852/(exp(-20.359146*pr - 7.4518167) + 1.0) + 0.00020245408;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end dynamicViscosity;

  redeclare function extends dBubbleEnthalpy_dPressure

  protected
    parameter AbsolutePressure pc=4901200;
    parameter SpecificEnthalpy h_ref=4e5;
    Real x;
  algorithm
    x:=sat.psat/pc;
    dhldp:=h_ref / pc * (4146.573/(exp(91.451248*x + 6.6011289) + 1.0) - 2292.3052/(exp(69.185605*x - 76.437727) + 1.0) + 1279.7798/(exp(21.801622 - 13.910316*x) + 1.0) + 1533.0444/(exp(5.6452119*x + 7.4171594) + 1.0) - 2386.0508/(exp(- 25.824527*x - 7.0233357) + 1.0) + 4678.6863);

  end dBubbleEnthalpy_dPressure;

  redeclare function extends dDewEnthalpy_dPressure

  protected
    parameter AbsolutePressure pc=4901200;
    parameter SpecificEnthalpy h_ref=4e5;
    Real x;
  algorithm
    x:=sat.psat/pc;
    dhvdp:=h_ref/pc*(1939.1795/(exp(70.776755*x - 78.381608) + 1.0) + 1515.5675/(exp(18.691841*x - 26.630322) + 1.0) + 1067.7974/(exp(2.958931*x - 10.920654) + 1.0) - 2191.2269/(exp(- 98.989675*x - 6.7533768) + 1.0) + 983.9616/(exp(6.6585351*x + 7.7408203) + 1.0) - 1404.8151/(exp(- 29.381977*x - 7.279896) + 1.0) - 926.47813);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end dDewEnthalpy_dPressure;

  redeclare function extends dBubbleDensity_dPressure
  protected
    Real x;
    parameter AbsolutePressure pc=4901200;
    parameter Density rhoc=459.03;
  algorithm
    x:=sat.psat/pc;
    ddldp:=rhoc / pc * (- 124.89469/(exp(31.9736/(exp(6.6066265*x - 9.2643555) + 1.0) + 39.394346/(exp(- 6.2008754*x - 2.8855421) + 1.0) - 65.927514) + 1.0) - 259.46317/(exp(188.4097/(exp(6.6066265*x - 9.2643555) + 1.0) + 232.76084/(exp(- 6.2008754*x - 2.8855421) + 1.0) - 405.67279) + 1.0) - 0.48403298);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end dBubbleDensity_dPressure;

  redeclare function extends dDewDensity_dPressure
  protected
    Real x;
    parameter AbsolutePressure pc=4901200;
    parameter Density rhoc=459.03;
  algorithm
    x:=sat.psat/pc;
    ddvdp:=rhoc / pc * (3005.254/(exp(25.244906 - 18.141441*x) + 1.0) - 4573.9243/(exp(69.348099*x - 76.176195) + 1.0) - 872.66696/(exp(3.5523568*x - 10.057065) + 1.0) + 5446.943);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end dDewDensity_dPressure;

  redeclare function extends bubbleEntropy
  protected
    Real p;
  algorithm
    p:=sat.psat;
    sl:=6537.3210414765195841877204677661/(exp(4.1504606546707407467090356968165 - 0.00000026882086611902591891035055252015*p) + 1.0) + 5711.179585869792626580760201629/(exp(-0.0000047733383332385529967848952663217*p - 3.3613509715459014546753411708511) + 1.0) + 5991.3325935223028752026663034312/(exp(-0.00000078145358720890553843915920837087*p -
      2.6860673677216475294051420354677) + 1.0) - 10612.747924970924558596977205757;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end bubbleEntropy;

  redeclare function extends dewEntropy
  protected
    Real p;
  algorithm
    p:=sat.psat;
    sv:=771.68234839099986027905426380517/(exp(0.00000064677916220004977374452634814104*p + 1.3882885534484957429347201630304) + 1.0) - 878.83859378768902892972283766367/(exp(-0.0000046749584597008177098477375356509*p - 2.0227154609816698767209559738363) + 1.0) - 825.70777469573380925400543788698/(exp(5.4245380859197011724041758533176 -
      0.00000078698759105517588386833372924335*p) + 1.0) + 2593.2814846909810717781343809111;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StartTime=620,
        StopTime=11830,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end dewEntropy;

  redeclare function extends specificEnthalpy_ps

  algorithm
    // Only vapor phase right now
    h:=2410647.6503256844342269983876052/(exp(0.0012943286865032919595990210301903*s - 0.000004556494230138293684619603827114*p - 6.4173498085461877494225051439403) + 1.0) + 836417.60096927973077908978305864/(exp(6.0028249492798649471728637911816 - 0.0020430331872715455857332103944182*s - 0.00000045589096001498377501485296745849*p) + 1.0) -
      1016718.62259152423299645461759/(exp(2.484581609829433121680398021209 - 0.00039138505216052310084171222435864*s - 0.00000042279126088122155374108011716905*p) + 1.0) + 1150134.5437367201695313749085421/(exp(0.0000008978457384509689446714681817584*p + 0.00073371828382245095292186069076974*s - 2.4733950842442032921886514097242) + 1.0) +
      1494863.9877968156823446130450643/(exp(2.6831096076965444224827701508903 - 0.00095422499225926154113150421017195*s - 0.00000081779609688898030808627873770531*p) + 1.0) - 3229848.0726115805134923093146146;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end specificEnthalpy_ps;

  redeclare function extends specificEntropy
  protected
    Real p=state.p;
    Real h=state.h;
  algorithm
    s:=109.6434865161011/(exp(0.00000124250852730057*p - 0.00001922566216648491*h + 7.677988828770944) + 1.0) - 33607.40325227791/(exp(2.333316588144453 - 0.0000002836436607334533*p - 0.000008911852400874945*h) + 1.0) - 53.8834457053386/(exp(0.0000297665771912091*h - 0.0000009769827799901408*p - 2.655447160182458) + 1.0) - 34783.37620958243/(exp(
      0.000008906994501832424*h + 0.0000002760352334200619*p - 2.30563072342852) + 1.0) - 758.765395187184/(exp(0.0000115572639463058*h - 0.0000001147687051024547*p - 4.654235635526474) + 1.0) + 354928.3415941631/(exp(0.000003343882343688367*p - 0.0000004708298756330904*h + 4.197410227478475) + 1.0) - 127868.9009697317/(exp(0.000003388236767065011*p -
      0.0000003630610378837811*h + 3.144859604011321) + 1.0) - 268.8287654232791/(exp(0.0000301093496974109*h - 0.0000002552545686062133*p - 1.982714631428223) + 1.0) + 35810.29172319622;

  end specificEntropy;

  redeclare function density_ph
    import DynamicVCC.Media.Utilities.phaseTransition;
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input Integer phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Density d "Density";
  protected
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Real x "Thermodynamic equilibrium quality";
    Density dl=bubbleDensity(sat);
    Density dv=dewDensity(sat);
    Density d_f "Liquid density model";
    Density d_g "Vapor denisty model";
    Density d_tp "Homogeneous model";
    Real pr=p/4901200 "Reduced pressure";
    Real hr=h/4e5 "Reduced enthalpy";
    Real w1;
    Real w2;
    Real w3;
  algorithm
    x:=(h - hl)/(hv - hl);
    (w1,w2,w3):=phaseTransition(x);

    d_f:=459.03*(2.6265285/(exp(5.9002435/(exp(0.043004485*pr - 3.2217628*hr + 1.1021958) + 1.0) - 46.408477/(exp(10.633616*hr - 1.0688275*pr - 12.17745) + 1.0) + 41.138912) + 1.0) - 0.80185018/(exp(458.41304/(exp(10.633616*hr - 1.0688275*pr - 12.17745) + 1.0) - 8.1497115/(exp(0.043004485*pr - 3.2217628*hr + 1.1021958) + 1.0) - 455.55522) + 1.0) + 1.1276562);

    d_g:=459.03*(0.57694216/(exp(4.2634059/(exp(2.0333087*pr - 1.3231126*hr + 0.78182062) + 1.0) - 255.8204/(exp(0.28641297*pr - 4.0813982*hr - 1.1775886) + 1.0) + 251.43574) + 1.0) + 3.2192717/(exp(224.11157/(exp(0.28641297*pr - 4.0813982*hr - 1.1775886) + 1.0) + 3.3953312/(exp(2.0333087*pr - 1.3231126*hr + 0.78182062) + 1.0) - 221.79698) + 1.0) - 0.45386158);

    d_tp:=1/(x/dv + (1 - x)/dl);

    d:=w1*d_f + w2*d_tp + w3*d_g;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end density_ph;

  redeclare function temperature_ph
    import DynamicVCC.Media.Utilities.phaseTransition;
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input Integer phase=0;
    output Temperature T "Temperature";
  protected
    parameter Temperature Tc=344.494;
    Temperature Tsat "Saturated temperature";
    SaturationProperties sat(psat=p,Tsat=T_default);
    Real Tl "Liquid temperature model";
    Real Tv "Vapor temperature model";
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Real x "Thermodynamic equilibrium quality";
    Real pr=p/4901200;
    Real hr=h/4e5;
    Real w1;
    Real w2;
    Real w3;
  algorithm
    Tsat:=saturationTemperature(p);

    x:=(h - hl)/(hv - hl);

    (w1,w2,w3):=phaseTransition(x);

    Tl:=Tc*(2.2485145/(exp(0.0077966862*pr - 1.5643748*hr + 0.18046657) + 1.0) - 35.83335/(exp(1.0617814*pr - 9.696756*hr + 13.716291) + 1.0) - 0.6549281);

    Tv:=Tc*(0.54723961/(exp(8.5240324*hr - 0.4401244*pr - 7.7964778) + 1.0) - 2.8093901/(exp(7.7218772*hr + 4.0511378*pr - 7.4533824) + 1.0) + 1.9631764/(exp(7.6406534*hr + 4.4883178*pr - 7.6190966) + 1.0) - 0.90562447/(exp(5.7546091*hr + 0.45799194*pr - 6.8938658) + 1.0) + 1.3873443);

    T:=w1*Tl + w2*Tsat + w3*Tv;



    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end temperature_ph;

  redeclare function extends density_derp_h
    import DynamicVCC.Media.Utilities.phaseTransition;
  protected
    parameter AbsolutePressure pc=4901200;
    parameter Density rhoc=459.03;
    parameter SpecificEnthalpy h_ref=4e5;
    Real pr=p/pc "Reduced pressure";
    Real hr "Reduced enthalpy";
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    Density rho=density(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Density rhof=bubbleDensity(sat);
    Density rhog=dewDensity(sat);
    Real xe;
    Real x "bounded in [0,1]";
    Real dxdp;
    DerDensityByPressure drhodp_l "subcooled liquid";
    DerDensityByPressure drhodp_v "superheated vapor";
    DerDensityByPressure drhodp_tp "two-phase";
    DerDensityByPressure drhodp_f=dBubbleDensity_dPressure(sat) "bubble line";
    DerDensityByPressure drhodp_g=dDewDensity_dPressure(sat) "dew line";
    DerEnthalpyByPressure dhdp_f=dBubbleEnthalpy_dPressure(sat) "bubble line";
    DerEnthalpyByPressure dhdp_g=dDewEnthalpy_dPressure(sat) "dew line";
    Real w1;
    Real w2;
    Real w3;
  algorithm
    xe:=(h - hl)/(hv - hl);
    (w1,w2,w3):=phaseTransition(xe);
    // bound x for two-phase model
    x:=(max(hl, min(hv, h)) - hl)/(hv - hl);
    dxdp:=(-x*dhdp_g - (1 - x)*dhdp_f)/(hv - hl);
    drhodp_tp:=x*(rho/rhog)^2*drhodp_g + (1 - x)*(rho/rhof)^2*drhodp_f - rho^2*(1/rhog - 1/rhof)*dxdp;
    // bound hr for subcooled liquid model
    hr:=min(hl, h)/h_ref;
    drhodp_l:=rhoc/pc*(5.2003707/(exp(1.0648072*pr - 14.931807*hr + 14.074691) + 1.0) - 0.80173175/(exp(66.370374*hr - 13.677169*pr - 43.993623) + 1.0) - 0.20988255/(exp(4.9608129*hr - 0.14107199*pr - 4.0730502) + 1.0) - 6.5015433/(exp(31.815383*hr - 1.384108*pr - 26.498629) + 1.0) - 9.9736239/(exp(0.8083154*pr - 30.12009*hr + 26.191205) + 1.0) + 7.5309909);
    // bound hr for superheated vapor model
    hr:=max(hv, h)/h_ref;
    drhodp_v:=rhoc/pc*(4.1059262/(exp(2.1189445*pr - 2.099238*hr + 0.50093889) + 1.0) + 4.4791415/(exp(2.4421226*hr - 1.988643*pr - 1.0285012) + 1.0) - 38.884364/(exp(1.2331259*pr - 15.357964*hr + 12.656595) + 1.0) + 15.758386/(exp(0.74654419*pr - 15.780767*hr + 13.612763) + 1.0) + 21.68091/(exp(1.7966177*pr - 15.473558*hr + 12.162137) + 1.0) - 2.4780512);
    // Fuzzy model
    ddph:=w1*drhodp_l + w2*drhodp_tp + w3*drhodp_v;




    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end density_derp_h;

  redeclare function extends density_derh_p
    import DynamicVCC.Media.Utilities.phaseTransition;
  protected
    parameter AbsolutePressure pc=4901200;
    parameter Density rhoc=459.03;
    parameter SpecificEnthalpy h_ref=4e5;
    Real pr=p/pc "Reduced pressure";
    Real hr "Reduced enthalpy";
    AbsolutePressure p=pressure(state);
    SpecificEnthalpy h=specificEnthalpy(state);
    Density rho=density(state);
    SaturationProperties sat(psat=p,Tsat=T_default);
    SpecificEnthalpy hl=bubbleEnthalpy(sat);
    SpecificEnthalpy hv=dewEnthalpy(sat);
    Density rhof=bubbleDensity(sat);
    Density rhog=dewDensity(sat);
    Real x;
    DerDensityByEnthalpy drhodh_l "Subcooled liquid";
    DerDensityByEnthalpy drhodh_v "Superheated vapor";
    DerDensityByEnthalpy drhodh_tp;
    Real w1;
    Real w2;
    Real w3;

  algorithm
    x:=(h - hl)/(hv - hl);
    (w1,w2,w3):=phaseTransition(x);
    // bound hr for subcooled liquid model
    hr:=min(hl, h)/h_ref;
    drhodh_l:=rhoc/h_ref*(69.164232/(exp(4.3616266*pr - 27.258951*hr + 19.212899) + 1.0) - 138.44746/(exp(4.0418744*pr - 26.494234*hr + 19.416731) + 1.0) - 0.99031126/(exp(0.35401042*pr - 5.1995607*hr + 2.7600425) + 1.0) - 71.305338/(exp(35.329515*hr - 3.8909459*pr - 28.148926) + 1.0) + 0.99205754/(exp(12.081155*hr - 0.8081261*pr - 7.5196917) + 1.0) + 68.566649);
    // bound hr for superheated vapor model
    hr:=max(hv, h)/h_ref;
    drhodh_v:=rhoc/h_ref*(376.38196/(exp(4.9606837*hr - 3.9219344*pr + 2.6950861) + 1.0) - 131.2105/(exp(6.2835703*hr - 1.437945*pr - 4.6712458) + 1.0) - 346.86334/(exp(4.5474352*hr - 4.6672042*pr + 4.3243156) + 1.0) + 151.02508/(exp(7.1540159*hr - 1.4149892*pr - 6.1902312) + 1.0) - 56.18676/(exp(7.931037*hr - 1.423085*pr - 7.350259) + 1.0) + 0.14847552);
    // two-phase model
    drhodh_tp:=-rho^2*(1/rhog - 1/rhof)/(hv - hl);
    // Fuzzy model
    ddhp:=min(-1e-8,w1*drhodh_l + w2*drhodh_tp + w3*drhodh_v);


    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end density_derh_p;

  redeclare function molarMass
    input ThermodynamicState state;
    output MolarMass MM;
  algorithm
    MM:=0.0725854;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end molarMass;

  redeclare function extends ThermalConductivity
    import DynamicVCC.Media.Utilities.phaseTransition;
  protected
    parameter AbsolutePressure pc=4901200;
    parameter SpecificEnthalpy h_ref=4e5;
    AbsolutePressure pr=p/pc;
    SpecificEnthalpy hr=h/h_ref;
  algorithm
    (w1,w2,w3):=phaseTransition(x);
    lambda_f:=2.7452654/(exp(33.297866*pr + 5.0305243) + 1.0) + 2.9112963/(exp(5.5810007*pr + 4.5187727) + 1.0) + 2.6969451/(exp(18.113615*pr - 23.907404) + 1.0) - 2.3016185/(exp(3.2420496 - 0.47271794*pr) + 1.0) - 2.513006;
    lambda_g:=3.1301689/(exp(32.248823 - 25.436557*pr) + 1.0) - 3.3395539/(exp(3.2363285*pr - 9.0380354) + 1.0) - 3.0109268/(exp(23.926491*pr + 7.0092493) + 1.0) + 3.1904209/(exp(-1.9819786*pr - 5.5774621) + 1.0) + 0.16968328;



    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
  end ThermalConductivity;
  annotation (Documentation(info="<html>
<p>Models of evaluating R410a properties are explicit in pressure and enthalpy (p,h), which are applicable within 3 to 45 bar and 106 to 472 kJ/kg.</p>
</html>"));
end R410a_NN;
