#include "specialized_solver/fast_monodomain_solver/tk2014.h"

#include <Vc/Vc>

Vc::double_v sqr(Vc::double_v x)
{
  return x*x;
}

Vc::double_v pow3(Vc::double_v x)
{
  return x*x*x;
}

Vc::double_v pow4(Vc::double_v x)
{
  return sqr(sqr(x));
}

double maxX = 0;
double minX = 0;

Vc::double_v exponential(Vc::double_v x)
{
  //return Vc::exp(x);
  // it was determined the x is always in the range [-12,+12]

  // exp(x) = lim n→∞ (1 + x/n)^n, we set n=1024
  x = 1.0 + x / 1024.;
  for (int i = 0; i < 10; i++)
  {
    x *= x;
  }
  return x;

  // relative error of this implementation:
  // x    rel error
  // 0    0
  // 1    0.00048784455634225593
  // 3    0.0043763626896140342
  // 5    0.012093715791500804
  // 9    0.038557535762274039
  // 12   0.067389808619653505
}

// new_slow_TK_2014_12_08
void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<57, 71, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
compute0DInstance(Vc::double_v states[], double currentTime, double timeStepWidth, bool stimulate)
{
  using Vc::double_v;

  // define constants
  const double constant0 = 0.58;
  const double constant1 = 2.79;
  const double constant2 = 150;
  const double constant3 = 0.000001;
  const double constant4 = 0.0025;
  const double constant5 = 0.0005;
  const double constant6 = 96485;
  const double constant7 = 559;
  const double constant8 = 559;
  const double constant9 = 0.00174;
  const double constant10 = 40229.885;
  const double constant11 = 40229.885;
  const double constant12 = 0.34;
  const double constant13 = -0.43;
  const double constant14 = 0.0081;
  const double constant15 = 0.288;
  const double constant16 = 0.0131;
  const double constant17 = 4.38;
  const double constant18 = 1.38;
  const double constant19 = 0.067;
  const double constant20 = -46;
  const double constant21 = -40;
  const double constant22 = -45;
  const double constant23 = 70;
  const double constant24 = -68;
  const double constant25 = -40;
  const double constant26 = 150;
  const double constant27 = 7.1;
  const double constant28 = 7.5;
  const double constant29 = 14.7;
  const double constant30 = 9;
  const double constant31 = 10;
  const double constant32 = 7;
  const double constant33 = 18;
  const double constant34 = 40;
  const double constant35 = 8314.41;
  const double constant36 = 293;
  const double constant37 = 3.275;
  const double constant38 = 10.8;
  const double constant39 = 134;
  const double constant40 = 1.85;
  const double constant41 = 0.4;
  const double constant42 = 950;
  const double constant43 = 1;
  const double constant44 = 1;
  const double constant45 = 13;
  const double constant46 = 10;
  const double constant47 = 0.0001656;
  const double constant48 = 70;
  const double constant49 = 0.1;
  const double constant50 = 1.0;
  const double constant51 = 0.45;
  const double constant52 = 0.1;
  const double constant53 = 0.1;
  const double constant54 = 0.0;
  const double constant55 = 0.002;
  const double constant56 = 1000;
  const double constant57 = 0.2;
  const double constant58 = 0.2;
  const double constant59 = 4.5;
  const double constant60 = -20;
  const double constant61 = 2.4375;
  const double constant62 = 1;
  const double constant63 = 0.00004;
  const double constant64 = 0.75;
  const double constant65 = 0.75;
  const double constant66 = 1.0;
  const double constant67 = 1.0;
  const double constant68 = 0.5;
  const double constant69 = 0.0885;
  const double constant70 = 0.115;
  const double constant71 = 140;
  const double constant72 = 0;
  const double constant73 = 0;
  const double constant74 = 1500;
  const double constant75 = 0;
  const double constant76 = 0;
  const double constant77 = 0.000004;
  const double constant78 = 0.005;
  const double constant79 = 31000;
  const double constant80 = 0.15;
  const double constant81 = 30;
  const double constant82 = 0.0015;
  const double constant83 = 0.15;
  const double constant84 = 0.375;
  const double constant85 = 1.5;
  const double constant86 = 0;
  const double constant87 = 0.15;
  const double constant88 = 0.15;
  const double constant89 = 0.05;
  const double constant90 = 0.5;
  const double constant91 = 5;
  const double constant92 = 0.08;
  const double constant93 = 0.06;
  const double constant94 = 0.04;
  const double constant95 = 0.00000394;
  const double constant96 = 0.00000362;
  const double constant97 = 1;
  const double constant98 = 0.0001;
  const double constant99 = 6;
  //const double constant100 = 0.05;  // this one is not used anywhere, hopefully this is not an error
  const double constant101 = 0.0;
  const double constant102 = 0.05;
  const double constant103 = 0.000107;
  const double constant104 = 0.0021;
  const double constant105 = 60;
  const double constant106 = 0.950000*constant66*3.14159265358979*constant68*constant68;
  const double constant107 = 0.0500000*constant66*3.14159265358979*constant68*constant68;
  const double constant108 = 0.0100000*constant106;
  const double constant109 = 0.990000*constant106;
  const double constant110 = 0.0100000*constant107;
  const double constant111 = 0.990000*constant107;

  // compute new rates, rhs(y_n)
  const double_v rate29 = (((( ( constant105*(states[18]+states[19]+states[20]+states[21]+states[22]))*((states[30] - states[29])/constant108) -  constant61*((states[29]/(states[29]+constant62))/constant108))+ constant63*((states[30] - states[29])/constant108))+ - constant64*((states[29] - states[31])/constant108))+- ( ( constant72*states[29])*((constant74+- states[34])+- states[36])+ - constant73*states[34]))+- ( ( constant80*states[29])*states[44]+ - constant81*states[40]);
  const double_v rate30 = ((( - ( constant105*(states[18]+states[19]+states[20]+states[21]+states[22]))*((states[30] - states[29])/constant110)+ constant61*((states[29]/(states[29]+constant62))/constant110))+ - constant63*((states[30] - states[29])/constant110))+ - constant65*((states[30] - states[32])/constant110))+- ( ( constant77*states[30])*(constant79 - states[38])+ - constant78*states[38]);
  const double_v rate32 = ((( constant61*((states[31]/(states[31]+constant62))/constant111)+ - constant63*((states[32]+- states[31])/constant111))+ constant65*((states[30]+- states[32])/constant111))+- ( ( constant77*states[32])*(constant79+- states[39])+ - constant78*states[39])) -  (1000.00/1.00000)*( constant97*( states[55]*(0.00100000/1.00000)*states[32] - constant99)*Vc::iif( states[55]*(0.00100000/1.00000)*states[32] - constant99>0.00000, double_v(Vc::One), double_v(Vc::Zero))*(0.00100000/1.00000)*states[55]*states[32] -  constant98*states[56]*(constant99 -  states[55]*(0.00100000/1.00000)*states[32])*Vc::iif(constant99 -  states[55]*(0.00100000/1.00000)*states[32]>0.00000, double_v(Vc::One), double_v(Vc::Zero)));
  const double_v rate34 =  ( constant72*states[29])*((constant74+- states[34])+- states[36])+ - constant73*states[34];
  const double_v rate35 =  ( constant72*states[31])*((constant74+- states[35])+- states[37])+ - constant73*states[35];
  const double_v rate36 =  ( constant75*(constant74+- states[34]+- states[36]))*states[46]+ - constant76*states[36];
  const double_v rate37 =  ( constant75*(constant74+- states[35]+- states[37]))*states[47]+ - constant76*states[37];
  const double_v rate38 =  ( constant77*states[30])*(constant79+- states[38])+ - constant78*states[38];
  const double_v rate39 =  ( constant77*states[32])*(constant79+- states[39])+ - constant78*states[39];
  const double_v rate40 = ( ( constant80*states[29])*states[44]+ - constant81*states[40])+ - constant84*((states[40]+- states[41])/constant108);
  const double_v rate41 = ( ( constant80*states[31])*states[45]+ - constant81*states[41])+ constant84*((states[40]+- states[41])/constant109);
  const double_v rate42 = ( ( constant82*states[46])*states[44]+ - constant83*states[42])+ - constant84*((states[42]+- states[43])/constant108);
  const double_v rate43 = ( ( constant82*states[47])*states[45]+ - constant83*states[43])+ constant84*((states[42]+- states[43])/constant109);
  const double_v rate44 = (- ( ( constant80*states[29])*states[44]+ - constant81*states[40])+- ( ( constant82*states[46])*states[44]+ - constant83*states[42]))+ - constant84*((states[44]+- states[45])/constant108);
  const double_v rate45 = (- ( ( constant80*states[31])*states[45]+ - constant81*states[41])+- ( ( constant82*states[47])*states[45]+ - constant83*states[43]))+ constant84*((states[44]+- states[45])/constant109);
  const double_v rate46 = (- ( ( constant75*(constant74+- states[34]+- states[36]))*states[46]+ - constant76*states[36])+- ( ( constant82*states[46])*states[44]+ - constant83*states[42]))+ - constant85*((states[46]+- states[47])/constant108);
  const double_v rate47 = (- ( ( constant75*(constant74+- states[35]+- states[37]))*states[47]+ - constant76*states[37])+- ( ( constant82*states[47])*states[45]+ - constant83*states[43]))+ constant85*((states[46]+- states[47])/constant109);
  const double_v rate48 = (( ( constant69*states[31])*states[33]+ - constant70*states[48])+ - constant88*states[48])+ constant89*states[51];
  const double_v rate50 = (((( constant69*states[31]*states[49]+ - constant70*states[50])+ constant86*states[33])+ - constant87*states[50])+ ( - constant69*states[31])*states[50])+ constant70*states[51];
  const double_v rate51 = ((((( constant69*states[31]*states[50]+ - constant70*states[51])+ constant88*states[48])+ - constant89*states[51])+ - constant90*states[51])+ constant91*states[52])+ constant94*states[53];
  const double_v rate52 = (( constant90*states[51]+ - constant91*states[52])+ constant93*states[53])+ - constant92*states[52];
  const double_v rate53 = ( - constant93*states[53]+ constant92*states[52])+ - constant94*states[53];
  const double_v rate54 =  (0.00100000/1.00000)*( constant92*states[52] -  constant93*states[53])+ -1.00000*constant95*states[54]+ -1.00000*constant96*((states[54] - states[55])/constant109);
  const double_v rate55 =  constant96*((states[54] - states[55])/constant111) -  1.00000*( constant97*( states[55]*(0.00100000/1.00000)*states[32] - constant99)*Vc::iif( states[55]*(0.00100000/1.00000)*states[32] - constant99>0.00000, double_v(Vc::One), double_v(Vc::Zero))*(0.00100000/1.00000)*states[55]*states[32] -  constant98*states[56]*(constant99 -  states[55]*(0.00100000/1.00000)*states[32])*Vc::iif(constant99 -  states[55]*(0.00100000/1.00000)*states[32]>0.00000, double_v(Vc::One), double_v(Vc::Zero)));
  const double_v rate56 =  1.00000*( constant97*( states[55]*(0.00100000/1.00000)*states[32] - constant99)*Vc::iif( states[55]*(0.00100000/1.00000)*states[32] - constant99>0.00000, double_v(Vc::One), double_v(Vc::Zero))*(0.00100000/1.00000)*states[55]*states[32] -  constant98*states[56]*(constant99 -  states[55]*(0.00100000/1.00000)*states[32])*Vc::iif(constant99 -  states[55]*(0.00100000/1.00000)*states[32]>0.00000, double_v(Vc::One), double_v(Vc::Zero)));
  const double_v algebraic13 = Vc::iif(constant71+- states[33]+- states[48]+- states[49]+- states[50]+- states[51]+- states[52]+- states[53]>0.00000 , constant71+- states[33]+- states[48]+- states[49]+- states[50]+- states[51]+- states[52]+- states[53] , double_v(Vc::Zero));
  const double_v rate31 = (((( - constant61*((states[31]/(states[31]+constant62))/constant109)+ constant63*((states[32]+- states[31])/constant109))+ constant64*((states[29] - states[31])/constant109))+- ((((((( constant69*states[31]*algebraic13+ - constant70*states[33])+ constant69*states[31]*states[33])+ - constant70*states[48])+ constant69*states[31]*states[49])+ - constant70*states[50])+ constant69*states[31]*states[50])+ - constant70*states[51]))+- ( ( constant72*states[31])*(constant74+- states[35]+- states[37])+ - constant73*states[35]))+- ( ( constant80*states[31])*states[45]+ - constant81*states[41]);
  const double_v rate33 = (((( ( constant69*states[31])*algebraic13+ - constant70*states[33])+ ( - constant69*states[31])*states[33])+ constant70*states[48])+ - constant86*states[33])+ constant87*states[50];
  const double_v rate49 = (( ( - constant69*states[31])*states[49]+ constant70*states[50])+ constant86*algebraic13)+ - constant87*states[49];
  const double_v algebraic12 =  ((( (states[52]/constant71)*constant101+ (states[53]/constant71)*constant102) - constant103)/constant104)*Vc::iif(constant67>=0.635000&&constant67<=0.850000,  (0.700000/(0.850000 - 0.635000))*(constant67 - 0.635000) , Vc::iif(constant67>0.850000&&constant67<=1.17000, 0.700000+ (0.300000/(1.17000 - 0.850000))*(constant67 - 0.850000), Vc::iif(constant67>1.17000&&constant67<=1.25500, 1.00000, Vc::iif(constant67>1.25500&&constant67<=1.97000, 1.00000 -  (1.00000/(1.97000 - 1.25500))*(constant67 - 1.25500), 0.00000))));
  const double_v rate28 = algebraic12;
  const double_v algebraic1 =  constant16*((states[0] - constant21)/(1.00000 - exponential(- ((states[0] - constant21)/constant32))));
  const double_v algebraic15 =  constant19*exponential(- ((states[0] - constant21)/constant34));
  const double_v rate8 =  algebraic1*(1.00000 - states[8]) -  algebraic15*states[8];
  const double_v algebraic2 = 1.00000/(1.00000+exponential((states[0] - constant25)/constant28));
  const double_v algebraic16 =  1000.00*exponential(- ((states[0]+40.0000)/25.7500));
  const double_v rate9 = (algebraic2 - states[9])/algebraic16;
  const double_v algebraic4 =  constant15*((states[0] - constant20)/(1.00000 - exponential(- ((states[0] - constant20)/constant31))));
  const double_v algebraic18 =  constant18*exponential(- ((states[0] - constant20)/constant33));
  const double_v rate10 =  algebraic4*(1.00000 - states[10]) -  algebraic18*states[10];
  const double_v algebraic3 =  constant14*exponential(- ((states[0] - constant22)/constant29));
  const double_v algebraic17 = constant17/(1.00000+exponential(- ((states[0] - constant22)/constant30)));
  const double_v rate11 =  algebraic3*(1.00000 - states[11]) -  algebraic17*states[11];
  const double_v algebraic5 = 1.00000/(1.00000+exponential((states[0] - constant24)/constant27));
  const double_v algebraic19 = 8571.00/(0.200000+ 5.65000*sqr((states[0]+constant48)/100.000));
  const double_v rate12 = (algebraic5 - states[12])/algebraic19;
  const double_v algebraic6 =  constant16*((states[1] - constant21)/(1.00000 - exponential(- ((states[1] - constant21)/constant32))));
  const double_v algebraic20 =  constant19*exponential(- ((states[1] - constant21)/constant34));
  const double_v rate13 =  algebraic6*(1.00000 - states[13]) -  algebraic20*states[13];
  const double_v algebraic7 = 1.00000/(1.00000+exponential((states[1] - constant25)/constant28));
  const double_v algebraic21 =  1.00000*exponential(- ((states[1]+40.0000)/25.7500));
  const double_v rate14 = (algebraic7 - states[14])/algebraic21;
  const double_v algebraic9 =  constant15*((states[1] - constant20)/(1.00000 - exponential(- ((states[1] - constant20)/constant31))));
  const double_v algebraic23 =  constant18*exponential(- ((states[1] - constant20)/constant33));
  const double_v rate15 =  algebraic9*(1.00000 - states[15]) -  algebraic23*states[15];
  const double_v algebraic8 =  constant14*exponential(- ((states[1] - constant22)/constant29));
  const double_v algebraic22 = constant17/(1.00000+exponential(- ((states[1] - constant22)/constant30)));
  const double_v rate16 =  algebraic8*(1.00000 - states[16]) -  algebraic22*states[16];
  const double_v algebraic10 = 1.00000/(1.00000+exponential((states[1] - constant24)/constant27));
  const double_v algebraic24 = 8571.00/(0.200000+ 5.65000*sqr((states[1]+constant48)/100.000));
  const double_v rate17 = (algebraic10 - states[17])/algebraic24;
  const double_v algebraic11 =  0.500000*constant58*exponential((states[1] - constant60)/( 8.00000*constant59));
  const double_v algebraic25 =  0.500000*constant58*exponential((constant60 - states[1])/( 8.00000*constant59));
  const double_v rate23 =  - constant55*states[23]+ constant56*states[18]+ -4.00000*algebraic11*states[23]+ algebraic25*states[24];
  const double_v rate18 =  constant55*states[23]+ - constant56*states[18]+( -4.00000*algebraic11*states[18])/constant57+ constant57*algebraic25*states[19];
  const double_v rate24 =  4.00000*algebraic11*states[23]+ - algebraic25*states[24]+( - constant55*states[24])/constant57+ constant57*constant56*states[19]+ -3.00000*algebraic11*states[24]+ 2.00000*algebraic25*states[25];
  const double_v rate19 = ( constant55*states[24])/constant57+ - constant56*constant57*states[19]+( 4.00000*algebraic11*states[18])/constant57+ - constant57*algebraic25*states[19]+( -3.00000*algebraic11*states[19])/constant57+ 2.00000*constant57*algebraic25*states[20];
  const double_v rate25 =  3.00000*algebraic11*states[24]+ -2.00000*algebraic25*states[25]+( - constant55*states[25])/pow(constant57, 2.00000)+ pow(constant57, 2.00000)*constant56*states[20]+ -2.00000*algebraic11*states[25]+ 3.00000*algebraic25*states[26];
  const double_v rate20 = ( 3.00000*algebraic11*states[19])/constant57+ -2.00000*constant57*algebraic25*states[20]+( constant55*states[25])/pow(constant57, 2.00000)+ - constant56*pow(constant57, 2.00000)*states[20]+( -2.00000*algebraic11*states[20])/constant57+ 3.00000*constant57*algebraic25*states[21];
  const double_v rate26 =  2.00000*algebraic11*states[25]+ -3.00000*algebraic25*states[26]+( - constant55*states[26])/pow(constant57, 3.00000)+ constant56*pow(constant57, 3.00000)*states[21]+ - algebraic11*states[26]+ 4.00000*algebraic25*states[27];
  const double_v rate21 = ( constant55*states[26])/pow(constant57, 3.00000)+ - constant56*pow(constant57, 3.00000)*states[21]+( 2.00000*algebraic11*states[20])/constant57+ -3.00000*algebraic25*constant57*states[21]+( - algebraic11*states[21])/constant57+ 4.00000*constant57*algebraic25*states[22];
  const double_v rate27 =  algebraic11*states[26]+ -4.00000*algebraic25*states[27]+( - constant55*states[27])/pow(constant57, 4.00000)+ constant56*pow(constant57, 4.00000)*states[22];
  const double_v rate22 = ( algebraic11*states[21])/constant57+ -4.00000*constant57*algebraic25*states[22]+( constant55*states[27])/pow(constant57, 4.00000)+ - constant56*pow(constant57, 4.00000)*states[22];
  const double_v algebraic31 =  states[0]*((states[3] -  states[4]*exponential(( -1.00000*constant6*states[0])/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*states[0])/( constant35*constant36))));
  const double_v algebraic14 =  (( constant35*constant36)/constant6)*log(states[4]/states[3]);
  const double_v algebraic37 =  states[4]*exponential( ( - constant41*algebraic14)*(constant6/( constant35*constant36)));
  const double_v algebraic38 =  constant40*(sqr(algebraic37)/(constant42+sqr(algebraic37)));
  const double_v algebraic39 = 1.00000 - 1.0/(1.00000+( constant43*(1.00000+sqr(algebraic37)/constant42))/( pow(constant46, 2.00000)*exponential(( 2.00000*(1.00000 - constant41)*states[0]*constant6)/( constant35*constant36))));
  const double_v algebraic40 =  algebraic38*algebraic39;
  const double_v algebraic41 =  algebraic40*Vc::iif(algebraic31>0.00000, double_v(Vc::One), double_v(Vc::Zero))*(algebraic31/50.0000);
  const double_v algebraic42 =  ( constant38*pow4(states[8]))*states[9];
  const double_v algebraic43 =  algebraic42*(algebraic31/50.0000);
  const double_v algebraic47 =  (1.00000/7.00000)*(exponential(states[7]/67.3000) - 1.00000);
  const double_v algebraic48 = 1.0/(1.00000+ 0.120000*exponential( -0.100000*states[0]*(constant6/( constant35*constant36)))+ 0.0400000*algebraic47*exponential(- ( states[0]*(constant6/( constant35*constant36)))));
  const double_v algebraic49 =  constant6*(constant47/( sqr(1.00000+constant44/states[4])*pow3(1.00000+constant45/states[5])));
  const double_v algebraic50 =  algebraic49*algebraic48;
  const double_v rate4 = (algebraic41+algebraic43+constant12+ - 2.00000*algebraic50)/( (1000.00/1.00000)*constant6*constant5)+(states[2] - states[4])/constant10;
  const double_v algebraic45 =  ( ( constant39*pow3(states[10]))*states[11])*states[12];
  const double_v algebraic44 =  states[0]*((states[5] -  states[7]*exponential(( -1.00000*constant6*states[0])/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*states[0])/( constant35*constant36))));
  const double_v algebraic46 =  algebraic45*(algebraic44/75.0000);
  const double_v rate7 = (algebraic46+constant13+ 3.00000*algebraic50)/( (1000.00/1.00000)*constant6*constant5)+(states[6] - states[7])/constant11;
  const double_v algebraic0 =  (1000.00/1.00000)*((states[0] - states[1])/constant2);
  const double_v algebraic27 = 156.500/(5.00000+exponential(( - constant6*algebraic14)/( constant35*constant36)));
  const double_v algebraic28 = 156.500 -  5.00000*algebraic27;
  const double_v algebraic34 =  states[0]*((algebraic27 -  algebraic28*exponential(( constant6*states[0])/( constant35*constant36)))/(1.00000 - exponential(( constant6*states[0])/( constant35*constant36))));
  const double_v algebraic33 = 1.00000/(1.00000+exponential((states[0] - constant23)/constant26));
  const double_v algebraic35 =  constant37*pow4(algebraic33);
  const double_v algebraic36 =  algebraic35*(algebraic34/45.0000);
  const double_v algebraic51 = algebraic36+algebraic41+algebraic43+algebraic46+algebraic50+- constant54;
  const double_v rate0 = - ((algebraic51+algebraic0)/constant0);
  const double_v algebraic32 =  states[1]*((states[3] -  states[2]*exponential(( -1.00000*constant6*states[1])/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*states[1])/( constant35*constant36))));
  const double_v algebraic26 =  (( constant35*constant36)/constant6)*log(states[2]/states[3]);
  const double_v algebraic56 =  states[2]*exponential( ( - constant41*algebraic26)*(constant6/( constant35*constant36)));
  const double_v algebraic57 =  constant40*(sqr(algebraic56)/(constant42+sqr(algebraic56)));
  const double_v algebraic58 = 1.00000 - 1.0/(1.00000+( constant43*(1.00000+sqr(algebraic56)/constant42))/( pow(constant46, 2.00000)*exponential(( 2.00000*(1.00000 - constant41)*states[1]*constant6)/( constant35*constant36))));
  const double_v algebraic59 =  algebraic57*algebraic58;
  const double_v algebraic60 =  constant50*algebraic59*(algebraic32/50.0000);
  const double_v algebraic61 =  ( constant38*pow4(states[13]))*states[14];
  const double_v algebraic62 =  constant51*algebraic61*(algebraic32/50.0000);
  const double_v algebraic66 =  (1.00000/7.00000)*(exponential(states[6]/67.3000) - 1.00000);
  const double_v algebraic67 = 1.0/(1.00000+ 0.120000*exponential( -0.100000*states[1]*(constant6/( constant35*constant36)))+ 0.0400000*algebraic66*exponential(- ( states[1]*(constant6/( constant35*constant36)))));
  const double_v algebraic68 =  constant6*(constant47/( sqr(1.00000+constant44/states[2])*pow3(1.00000+constant45/states[5])));
  const double_v algebraic69 =  constant53*algebraic68*algebraic67;
  const double_v rate3 =  - constant9*((algebraic60+algebraic62+constant12+ - 2.00000*algebraic69)/( (1000.00/1.00000)*constant6*constant3)) - (algebraic41+algebraic43+constant12+ -2.00000*algebraic50)/( (1000.00/1.00000)*constant6*constant4);
  const double_v rate2 = (algebraic60+algebraic62+constant12+ - 2.00000*algebraic69)/( (1000.00/1.00000)*constant6*constant3) - (states[2] - states[4])/constant7;
  const double_v algebraic64 =  ( ( constant39*pow3(states[15]))*states[16])*states[17];
  const double_v algebraic63 =  states[1]*((states[5] -  states[6]*exponential(( -1.00000*constant6*states[1])/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*states[1])/( constant35*constant36))));
  const double_v algebraic65 =  constant52*algebraic64*(algebraic63/75.0000);
  const double_v rate5 =  - constant9*((algebraic65+constant13+ 3.00000*algebraic69)/( (1000.00/1.00000)*constant6*constant3)) - (algebraic46+constant13+ 3.00000*algebraic50)/( (1000.00/1.00000)*constant6*constant4);
  const double_v rate6 = (algebraic65+constant13+ 3.00000*algebraic69)/( (1000.00/1.00000)*constant6*constant3) - (states[6] - states[7])/constant8;
  const double_v algebraic29 = 156.500/(5.00000+exponential(( - constant6*algebraic26)/( constant35*constant36)));
  const double_v algebraic30 = 156.500 -  5.00000*algebraic29;
  const double_v algebraic53 =  states[1]*((algebraic29 -  algebraic30*exponential(( constant6*states[1])/( constant35*constant36)))/(1.00000 - exponential(( constant6*states[1])/( constant35*constant36))));
  const double_v algebraic52 = 1.00000/(1.00000+exponential((states[1] - constant23)/constant26));
  const double_v algebraic54 =  constant37*pow4(algebraic52);
  const double_v algebraic55 =  constant49*algebraic54*(algebraic53/45.0000);
  const double_v algebraic70 = algebraic55+algebraic60+algebraic62+algebraic65+algebraic69;
  const double_v rate1 = - ((algebraic70 - algebraic0/constant1)/constant0);

  // intermediate step
  // compute y* = y_n + dt*rhs(y_n), y_n = state, rhs(y_n) = rate, y* = intermediateState
  double_v intermediateState0 = states[0] + timeStepWidth*rate0;
  const double_v intermediateState1 = states[1] + timeStepWidth*rate1;
  const double_v intermediateState2 = states[2] + timeStepWidth*rate2;
  const double_v intermediateState3 = states[3] + timeStepWidth*rate3;
  const double_v intermediateState4 = states[4] + timeStepWidth*rate4;
  const double_v intermediateState5 = states[5] + timeStepWidth*rate5;
  const double_v intermediateState6 = states[6] + timeStepWidth*rate6;
  const double_v intermediateState7 = states[7] + timeStepWidth*rate7;
  const double_v intermediateState8 = states[8] + timeStepWidth*rate8;
  const double_v intermediateState9 = states[9] + timeStepWidth*rate9;
  const double_v intermediateState10 = states[10] + timeStepWidth*rate10;
  const double_v intermediateState11 = states[11] + timeStepWidth*rate11;
  const double_v intermediateState12 = states[12] + timeStepWidth*rate12;
  const double_v intermediateState13 = states[13] + timeStepWidth*rate13;
  const double_v intermediateState14 = states[14] + timeStepWidth*rate14;
  const double_v intermediateState15 = states[15] + timeStepWidth*rate15;
  const double_v intermediateState16 = states[16] + timeStepWidth*rate16;
  const double_v intermediateState17 = states[17] + timeStepWidth*rate17;
  const double_v intermediateState18 = states[18] + timeStepWidth*rate18;
  const double_v intermediateState19 = states[19] + timeStepWidth*rate19;
  const double_v intermediateState20 = states[20] + timeStepWidth*rate20;
  const double_v intermediateState21 = states[21] + timeStepWidth*rate21;
  const double_v intermediateState22 = states[22] + timeStepWidth*rate22;
  const double_v intermediateState23 = states[23] + timeStepWidth*rate23;
  const double_v intermediateState24 = states[24] + timeStepWidth*rate24;
  const double_v intermediateState25 = states[25] + timeStepWidth*rate25;
  const double_v intermediateState26 = states[26] + timeStepWidth*rate26;
  const double_v intermediateState27 = states[27] + timeStepWidth*rate27;
//  const double_v intermediateState28 = states[28] + timeStepWidth*rate28;
  const double_v intermediateState29 = states[29] + timeStepWidth*rate29;
  const double_v intermediateState30 = states[30] + timeStepWidth*rate30;
  const double_v intermediateState31 = states[31] + timeStepWidth*rate31;
  const double_v intermediateState32 = states[32] + timeStepWidth*rate32;
  const double_v intermediateState33 = states[33] + timeStepWidth*rate33;
  const double_v intermediateState34 = states[34] + timeStepWidth*rate34;
  const double_v intermediateState35 = states[35] + timeStepWidth*rate35;
  const double_v intermediateState36 = states[36] + timeStepWidth*rate36;
  const double_v intermediateState37 = states[37] + timeStepWidth*rate37;
  const double_v intermediateState38 = states[38] + timeStepWidth*rate38;
  const double_v intermediateState39 = states[39] + timeStepWidth*rate39;
  const double_v intermediateState40 = states[40] + timeStepWidth*rate40;
  const double_v intermediateState41 = states[41] + timeStepWidth*rate41;
  const double_v intermediateState42 = states[42] + timeStepWidth*rate42;
  const double_v intermediateState43 = states[43] + timeStepWidth*rate43;
  const double_v intermediateState44 = states[44] + timeStepWidth*rate44;
  const double_v intermediateState45 = states[45] + timeStepWidth*rate45;
  const double_v intermediateState46 = states[46] + timeStepWidth*rate46;
  const double_v intermediateState47 = states[47] + timeStepWidth*rate47;
  const double_v intermediateState48 = states[48] + timeStepWidth*rate48;
  const double_v intermediateState49 = states[49] + timeStepWidth*rate49;
  const double_v intermediateState50 = states[50] + timeStepWidth*rate50;
  const double_v intermediateState51 = states[51] + timeStepWidth*rate51;
  const double_v intermediateState52 = states[52] + timeStepWidth*rate52;
  const double_v intermediateState53 = states[53] + timeStepWidth*rate53;
  const double_v intermediateState54 = states[54] + timeStepWidth*rate54;
  const double_v intermediateState55 = states[55] + timeStepWidth*rate55;
  const double_v intermediateState56 = states[56] + timeStepWidth*rate56;


  // if stimulation, set value of Vm (state0)
  if (stimulate)
  {
    for (int i = 0; i < std::min(3,(int)Vc::double_v::Size); i++)
    {
      intermediateState0[i] = 20.0;
    }
  }

  // compute new rates, rhs(y*)
  const double_v intermediateRate29 = (((( ( constant105*(intermediateState18+intermediateState19+intermediateState20+intermediateState21+intermediateState22))*((intermediateState30 - intermediateState29)/constant108) -  constant61*((intermediateState29/(intermediateState29+constant62))/constant108))+ constant63*((intermediateState30 - intermediateState29)/constant108))+ - constant64*((intermediateState29 - intermediateState31)/constant108))+- ( ( constant72*intermediateState29)*((constant74+- intermediateState34)+- intermediateState36)+ - constant73*intermediateState34))+- ( ( constant80*intermediateState29)*intermediateState44+ - constant81*intermediateState40);
  const double_v intermediateRate30 = ((( - ( constant105*(intermediateState18+intermediateState19+intermediateState20+intermediateState21+intermediateState22))*((intermediateState30 - intermediateState29)/constant110)+ constant61*((intermediateState29/(intermediateState29+constant62))/constant110))+ - constant63*((intermediateState30 - intermediateState29)/constant110))+ - constant65*((intermediateState30 - intermediateState32)/constant110))+- ( ( constant77*intermediateState30)*(constant79 - intermediateState38)+ - constant78*intermediateState38);
  const double_v intermediateRate32 = ((( constant61*((intermediateState31/(intermediateState31+constant62))/constant111)+ - constant63*((intermediateState32+- intermediateState31)/constant111))+ constant65*((intermediateState30+- intermediateState32)/constant111))+- ( ( constant77*intermediateState32)*(constant79+- intermediateState39)+ - constant78*intermediateState39)) -  (1000.00/1.00000)*( constant97*( intermediateState55*(0.00100000/1.00000)*intermediateState32 - constant99)*Vc::iif( intermediateState55*(0.00100000/1.00000)*intermediateState32 - constant99>0.00000, double_v(Vc::One), double_v(Vc::Zero))*(0.00100000/1.00000)*intermediateState55*intermediateState32 -  constant98*intermediateState56*(constant99 -  intermediateState55*(0.00100000/1.00000)*intermediateState32)*Vc::iif(constant99 -  intermediateState55*(0.00100000/1.00000)*intermediateState32>0.00000, double_v(Vc::One), double_v(Vc::Zero)));
  const double_v intermediateRate34 =  ( constant72*intermediateState29)*((constant74+- intermediateState34)+- intermediateState36)+ - constant73*intermediateState34;
  const double_v intermediateRate35 =  ( constant72*intermediateState31)*((constant74+- intermediateState35)+- intermediateState37)+ - constant73*intermediateState35;
  const double_v intermediateRate36 =  ( constant75*(constant74+- intermediateState34+- intermediateState36))*intermediateState46+ - constant76*intermediateState36;
  const double_v intermediateRate37 =  ( constant75*(constant74+- intermediateState35+- intermediateState37))*intermediateState47+ - constant76*intermediateState37;
  const double_v intermediateRate38 =  ( constant77*intermediateState30)*(constant79+- intermediateState38)+ - constant78*intermediateState38;
  const double_v intermediateRate39 =  ( constant77*intermediateState32)*(constant79+- intermediateState39)+ - constant78*intermediateState39;
  const double_v intermediateRate40 = ( ( constant80*intermediateState29)*intermediateState44+ - constant81*intermediateState40)+ - constant84*((intermediateState40+- intermediateState41)/constant108);
  const double_v intermediateRate41 = ( ( constant80*intermediateState31)*intermediateState45+ - constant81*intermediateState41)+ constant84*((intermediateState40+- intermediateState41)/constant109);
  const double_v intermediateRate42 = ( ( constant82*intermediateState46)*intermediateState44+ - constant83*intermediateState42)+ - constant84*((intermediateState42+- intermediateState43)/constant108);
  const double_v intermediateRate43 = ( ( constant82*intermediateState47)*intermediateState45+ - constant83*intermediateState43)+ constant84*((intermediateState42+- intermediateState43)/constant109);
  const double_v intermediateRate44 = (- ( ( constant80*intermediateState29)*intermediateState44+ - constant81*intermediateState40)+- ( ( constant82*intermediateState46)*intermediateState44+ - constant83*intermediateState42))+ - constant84*((intermediateState44+- intermediateState45)/constant108);
  const double_v intermediateRate45 = (- ( ( constant80*intermediateState31)*intermediateState45+ - constant81*intermediateState41)+- ( ( constant82*intermediateState47)*intermediateState45+ - constant83*intermediateState43))+ constant84*((intermediateState44+- intermediateState45)/constant109);
  const double_v intermediateRate46 = (- ( ( constant75*(constant74+- intermediateState34+- intermediateState36))*intermediateState46+ - constant76*intermediateState36)+- ( ( constant82*intermediateState46)*intermediateState44+ - constant83*intermediateState42))+ - constant85*((intermediateState46+- intermediateState47)/constant108);
  const double_v intermediateRate47 = (- ( ( constant75*(constant74+- intermediateState35+- intermediateState37))*intermediateState47+ - constant76*intermediateState37)+- ( ( constant82*intermediateState47)*intermediateState45+ - constant83*intermediateState43))+ constant85*((intermediateState46+- intermediateState47)/constant109);
  const double_v intermediateRate48 = (( ( constant69*intermediateState31)*intermediateState33+ - constant70*intermediateState48)+ - constant88*intermediateState48)+ constant89*intermediateState51;
  const double_v intermediateRate50 = (((( constant69*intermediateState31*intermediateState49+ - constant70*intermediateState50)+ constant86*intermediateState33)+ - constant87*intermediateState50)+ ( - constant69*intermediateState31)*intermediateState50)+ constant70*intermediateState51;
  const double_v intermediateRate51 = ((((( constant69*intermediateState31*intermediateState50+ - constant70*intermediateState51)+ constant88*intermediateState48)+ - constant89*intermediateState51)+ - constant90*intermediateState51)+ constant91*intermediateState52)+ constant94*intermediateState53;
  const double_v intermediateRate52 = (( constant90*intermediateState51+ - constant91*intermediateState52)+ constant93*intermediateState53)+ - constant92*intermediateState52;
  const double_v intermediateRate53 = ( - constant93*intermediateState53+ constant92*intermediateState52)+ - constant94*intermediateState53;
  const double_v intermediateRate54 =  (0.00100000/1.00000)*( constant92*intermediateState52 -  constant93*intermediateState53)+ -1.00000*constant95*intermediateState54+ -1.00000*constant96*((intermediateState54 - intermediateState55)/constant109);
  const double_v intermediateRate55 =  constant96*((intermediateState54 - intermediateState55)/constant111) -  1.00000*( constant97*( intermediateState55*(0.00100000/1.00000)*intermediateState32 - constant99)*Vc::iif( intermediateState55*(0.00100000/1.00000)*intermediateState32 - constant99>0.00000, double_v(Vc::One), double_v(Vc::Zero))*(0.00100000/1.00000)*intermediateState55*intermediateState32 -  constant98*intermediateState56*(constant99 -  intermediateState55*(0.00100000/1.00000)*intermediateState32)*Vc::iif(constant99 -  intermediateState55*(0.00100000/1.00000)*intermediateState32>0.00000, double_v(Vc::One), double_v(Vc::Zero)));
  const double_v intermediateRate56 =  1.00000*( constant97*( intermediateState55*(0.00100000/1.00000)*intermediateState32 - constant99)*Vc::iif( intermediateState55*(0.00100000/1.00000)*intermediateState32 - constant99>0.00000, double_v(Vc::One) , double_v(Vc::Zero))*(0.00100000/1.00000)*intermediateState55*intermediateState32 -  constant98*intermediateState56*(constant99 -  intermediateState55*(0.00100000/1.00000)*intermediateState32)*Vc::iif(constant99 -  intermediateState55*(0.00100000/1.00000)*intermediateState32>0.00000 , double_v(Vc::One) , double_v(Vc::Zero)));
  const double_v intermediateAlgebraic13 = Vc::iif(constant71+- intermediateState33+- intermediateState48+- intermediateState49+- intermediateState50+- intermediateState51+- intermediateState52+- intermediateState53>0.00000 , constant71+- intermediateState33+- intermediateState48+- intermediateState49+- intermediateState50+- intermediateState51+- intermediateState52+- intermediateState53 , double_v(Vc::Zero));
  const double_v intermediateRate31 = (((( - constant61*((intermediateState31/(intermediateState31+constant62))/constant109)+ constant63*((intermediateState32+- intermediateState31)/constant109))+ constant64*((intermediateState29 - intermediateState31)/constant109))+- ((((((( constant69*intermediateState31*intermediateAlgebraic13+ - constant70*intermediateState33)+ constant69*intermediateState31*intermediateState33)+ - constant70*intermediateState48)+ constant69*intermediateState31*intermediateState49)+ - constant70*intermediateState50)+ constant69*intermediateState31*intermediateState50)+ - constant70*intermediateState51))+- ( ( constant72*intermediateState31)*(constant74+- intermediateState35+- intermediateState37)+ - constant73*intermediateState35))+- ( ( constant80*intermediateState31)*intermediateState45+ - constant81*intermediateState41);
  const double_v intermediateRate33 = (((( ( constant69*intermediateState31)*intermediateAlgebraic13+ - constant70*intermediateState33)+ ( - constant69*intermediateState31)*intermediateState33)+ constant70*intermediateState48)+ - constant86*intermediateState33)+ constant87*intermediateState50;
  const double_v intermediateRate49 = (( ( - constant69*intermediateState31)*intermediateState49+ constant70*intermediateState50)+ constant86*intermediateAlgebraic13)+ - constant87*intermediateState49;
  const double_v intermediateAlgebraic12 =  ((( (intermediateState52/constant71)*constant101+ (intermediateState53/constant71)*constant102) - constant103)/constant104)*Vc::iif(constant67>=0.635000&&constant67<=0.850000 ,  (0.700000/(0.850000 - 0.635000))*(constant67 - 0.635000) , Vc::iif(constant67>0.850000&&constant67<=1.17000 , 0.700000+ (0.300000/(1.17000 - 0.850000))*(constant67 - 0.850000) , Vc::iif(constant67>1.17000&&constant67<=1.25500 , 1.00000 , Vc::iif(constant67>1.25500&&constant67<=1.97000 , 1.00000 -  (1.00000/(1.97000 - 1.25500))*(constant67 - 1.25500) , 0.00000))));
  const double_v intermediateRate28 = intermediateAlgebraic12;
  const double_v intermediateAlgebraic1 =  constant16*((intermediateState0 - constant21)/(1.00000 - exponential(- ((intermediateState0 - constant21)/constant32))));
  const double_v intermediateAlgebraic15 =  constant19*exponential(- ((intermediateState0 - constant21)/constant34));
  const double_v intermediateRate8 =  intermediateAlgebraic1*(1.00000 - intermediateState8) -  intermediateAlgebraic15*intermediateState8;
  const double_v intermediateAlgebraic2 = 1.00000/(1.00000+exponential((intermediateState0 - constant25)/constant28));
  const double_v intermediateAlgebraic16 =  1000.00*exponential(- ((intermediateState0+40.0000)/25.7500));
  const double_v intermediateRate9 = (intermediateAlgebraic2 - intermediateState9)/intermediateAlgebraic16;
  const double_v intermediateAlgebraic4 =  constant15*((intermediateState0 - constant20)/(1.00000 - exponential(- ((intermediateState0 - constant20)/constant31))));
  const double_v intermediateAlgebraic18 =  constant18*exponential(- ((intermediateState0 - constant20)/constant33));
  const double_v intermediateRate10 =  intermediateAlgebraic4*(1.00000 - intermediateState10) -  intermediateAlgebraic18*intermediateState10;
  const double_v intermediateAlgebraic3 =  constant14*exponential(- ((intermediateState0 - constant22)/constant29));
  const double_v intermediateAlgebraic17 = constant17/(1.00000+exponential(- ((intermediateState0 - constant22)/constant30)));
  const double_v intermediateRate11 =  intermediateAlgebraic3*(1.00000 - intermediateState11) -  intermediateAlgebraic17*intermediateState11;
  const double_v intermediateAlgebraic5 = 1.00000/(1.00000+exponential((intermediateState0 - constant24)/constant27));
  const double_v intermediateAlgebraic19 = 8571.00/(0.200000+ 5.65000*sqr((intermediateState0+constant48)/100.000));
  const double_v intermediateRate12 = (intermediateAlgebraic5 - intermediateState12)/intermediateAlgebraic19;
  const double_v intermediateAlgebraic6 =  constant16*((intermediateState1 - constant21)/(1.00000 - exponential(- ((intermediateState1 - constant21)/constant32))));
  const double_v intermediateAlgebraic20 =  constant19*exponential(- ((intermediateState1 - constant21)/constant34));
  const double_v intermediateRate13 =  intermediateAlgebraic6*(1.00000 - intermediateState13) -  intermediateAlgebraic20*intermediateState13;
  const double_v intermediateAlgebraic7 = 1.00000/(1.00000+exponential((intermediateState1 - constant25)/constant28));
  const double_v intermediateAlgebraic21 =  1.00000*exponential(- ((intermediateState1+40.0000)/25.7500));
  const double_v intermediateRate14 = (intermediateAlgebraic7 - intermediateState14)/intermediateAlgebraic21;
  const double_v intermediateAlgebraic9 =  constant15*((intermediateState1 - constant20)/(1.00000 - exponential(- ((intermediateState1 - constant20)/constant31))));
  const double_v intermediateAlgebraic23 =  constant18*exponential(- ((intermediateState1 - constant20)/constant33));
  const double_v intermediateRate15 =  intermediateAlgebraic9*(1.00000 - intermediateState15) -  intermediateAlgebraic23*intermediateState15;
  const double_v intermediateAlgebraic8 =  constant14*exponential(- ((intermediateState1 - constant22)/constant29));
  const double_v intermediateAlgebraic22 = constant17/(1.00000+exponential(- ((intermediateState1 - constant22)/constant30)));
  const double_v intermediateRate16 =  intermediateAlgebraic8*(1.00000 - intermediateState16) -  intermediateAlgebraic22*intermediateState16;
  const double_v intermediateAlgebraic10 = 1.00000/(1.00000+exponential((intermediateState1 - constant24)/constant27));
  const double_v intermediateAlgebraic24 = 8571.00/(0.200000+ 5.65000*sqr((intermediateState1+constant48)/100.000));
  const double_v intermediateRate17 = (intermediateAlgebraic10 - intermediateState17)/intermediateAlgebraic24;
  const double_v intermediateAlgebraic11 =  0.500000*constant58*exponential((intermediateState1 - constant60)/( 8.00000*constant59));
  const double_v intermediateAlgebraic25 =  0.500000*constant58*exponential((constant60 - intermediateState1)/( 8.00000*constant59));
  const double_v intermediateRate23 =  - constant55*intermediateState23+ constant56*intermediateState18+ -4.00000*intermediateAlgebraic11*intermediateState23+ intermediateAlgebraic25*intermediateState24;
  const double_v intermediateRate18 =  constant55*intermediateState23+ - constant56*intermediateState18+( -4.00000*intermediateAlgebraic11*intermediateState18)/constant57+ constant57*intermediateAlgebraic25*intermediateState19;
  const double_v intermediateRate24 =  4.00000*intermediateAlgebraic11*intermediateState23+ - intermediateAlgebraic25*intermediateState24+( - constant55*intermediateState24)/constant57+ constant57*constant56*intermediateState19+ -3.00000*intermediateAlgebraic11*intermediateState24+ 2.00000*intermediateAlgebraic25*intermediateState25;
  const double_v intermediateRate19 = ( constant55*intermediateState24)/constant57+ - constant56*constant57*intermediateState19+( 4.00000*intermediateAlgebraic11*intermediateState18)/constant57+ - constant57*intermediateAlgebraic25*intermediateState19+( -3.00000*intermediateAlgebraic11*intermediateState19)/constant57+ 2.00000*constant57*intermediateAlgebraic25*intermediateState20;
  const double_v intermediateRate25 =  3.00000*intermediateAlgebraic11*intermediateState24+ -2.00000*intermediateAlgebraic25*intermediateState25+( - constant55*intermediateState25)/pow(constant57, 2.00000)+ pow(constant57, 2.00000)*constant56*intermediateState20+ -2.00000*intermediateAlgebraic11*intermediateState25+ 3.00000*intermediateAlgebraic25*intermediateState26;
  const double_v intermediateRate20 = ( 3.00000*intermediateAlgebraic11*intermediateState19)/constant57+ -2.00000*constant57*intermediateAlgebraic25*intermediateState20+( constant55*intermediateState25)/pow(constant57, 2.00000)+ - constant56*pow(constant57, 2.00000)*intermediateState20+( -2.00000*intermediateAlgebraic11*intermediateState20)/constant57+ 3.00000*constant57*intermediateAlgebraic25*intermediateState21;
  const double_v intermediateRate26 =  2.00000*intermediateAlgebraic11*intermediateState25+ -3.00000*intermediateAlgebraic25*intermediateState26+( - constant55*intermediateState26)/pow(constant57, 3.00000)+ constant56*pow(constant57, 3.00000)*intermediateState21+ - intermediateAlgebraic11*intermediateState26+ 4.00000*intermediateAlgebraic25*intermediateState27;
  const double_v intermediateRate21 = ( constant55*intermediateState26)/pow(constant57, 3.00000)+ - constant56*pow(constant57, 3.00000)*intermediateState21+( 2.00000*intermediateAlgebraic11*intermediateState20)/constant57+ -3.00000*intermediateAlgebraic25*constant57*intermediateState21+( - intermediateAlgebraic11*intermediateState21)/constant57+ 4.00000*constant57*intermediateAlgebraic25*intermediateState22;
  const double_v intermediateRate27 =  intermediateAlgebraic11*intermediateState26+ -4.00000*intermediateAlgebraic25*intermediateState27+( - constant55*intermediateState27)/pow(constant57, 4.00000)+ constant56*pow(constant57, 4.00000)*intermediateState22;
  const double_v intermediateRate22 = ( intermediateAlgebraic11*intermediateState21)/constant57+ -4.00000*constant57*intermediateAlgebraic25*intermediateState22+( constant55*intermediateState27)/pow(constant57, 4.00000)+ - constant56*pow(constant57, 4.00000)*intermediateState22;
  const double_v intermediateAlgebraic31 =  intermediateState0*((intermediateState3 -  intermediateState4*exponential(( -1.00000*constant6*intermediateState0)/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*intermediateState0)/( constant35*constant36))));
  const double_v intermediateAlgebraic14 =  (( constant35*constant36)/constant6)*log(intermediateState4/intermediateState3);
  const double_v intermediateAlgebraic37 =  intermediateState4*exponential( ( - constant41*intermediateAlgebraic14)*(constant6/( constant35*constant36)));
  const double_v intermediateAlgebraic38 =  constant40*(sqr(intermediateAlgebraic37)/(constant42+sqr(intermediateAlgebraic37)));
  const double_v intermediateAlgebraic39 = 1.00000 - 1.0/(1.00000+( constant43*(1.00000+sqr(intermediateAlgebraic37)/constant42))/( pow(constant46, 2.00000)*exponential(( 2.00000*(1.00000 - constant41)*intermediateState0*constant6)/( constant35*constant36))));
  const double_v intermediateAlgebraic40 =  intermediateAlgebraic38*intermediateAlgebraic39;
  const double_v intermediateAlgebraic41 =  intermediateAlgebraic40*Vc::iif(intermediateAlgebraic31>0.00000 , double_v(Vc::One) , double_v(Vc::Zero))*(intermediateAlgebraic31/50.0000);
  const double_v intermediateAlgebraic42 =  ( constant38*pow4(intermediateState8))*intermediateState9;
  const double_v intermediateAlgebraic43 =  intermediateAlgebraic42*(intermediateAlgebraic31/50.0000);
  const double_v intermediateAlgebraic47 =  (1.00000/7.00000)*(exponential(intermediateState7/67.3000) - 1.00000);
  const double_v intermediateAlgebraic48 = 1.0/(1.00000+ 0.120000*exponential( -0.100000*intermediateState0*(constant6/( constant35*constant36)))+ 0.0400000*intermediateAlgebraic47*exponential(- ( intermediateState0*(constant6/( constant35*constant36)))));
  const double_v intermediateAlgebraic49 =  constant6*(constant47/( sqr(1.00000+constant44/intermediateState4)*pow3(1.00000+constant45/intermediateState5)));
  const double_v intermediateAlgebraic50 =  intermediateAlgebraic49*intermediateAlgebraic48;
  const double_v intermediateRate4 = (intermediateAlgebraic41+intermediateAlgebraic43+constant12+ - 2.00000*intermediateAlgebraic50)/( (1000.00/1.00000)*constant6*constant5)+(intermediateState2 - intermediateState4)/constant10;
  const double_v intermediateAlgebraic45 =  ( ( constant39*pow3(intermediateState10))*intermediateState11)*intermediateState12;
  const double_v intermediateAlgebraic44 =  intermediateState0*((intermediateState5 -  intermediateState7*exponential(( -1.00000*constant6*intermediateState0)/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*intermediateState0)/( constant35*constant36))));
  const double_v intermediateAlgebraic46 =  intermediateAlgebraic45*(intermediateAlgebraic44/75.0000);
  const double_v intermediateRate7 = (intermediateAlgebraic46+constant13+ 3.00000*intermediateAlgebraic50)/( (1000.00/1.00000)*constant6*constant5)+(intermediateState6 - intermediateState7)/constant11;
  const double_v intermediateAlgebraic0 =  (1000.00/1.00000)*((intermediateState0 - intermediateState1)/constant2);
  const double_v intermediateAlgebraic27 = 156.500/(5.00000+exponential(( - constant6*intermediateAlgebraic14)/( constant35*constant36)));
  const double_v intermediateAlgebraic28 = 156.500 -  5.00000*intermediateAlgebraic27;
  const double_v intermediateAlgebraic34 =  intermediateState0*((intermediateAlgebraic27 -  intermediateAlgebraic28*exponential(( constant6*intermediateState0)/( constant35*constant36)))/(1.00000 - exponential(( constant6*intermediateState0)/( constant35*constant36))));
  const double_v intermediateAlgebraic33 = 1.00000/(1.00000+exponential((intermediateState0 - constant23)/constant26));
  const double_v intermediateAlgebraic35 =  constant37*pow4(intermediateAlgebraic33);
  const double_v intermediateAlgebraic36 =  intermediateAlgebraic35*(intermediateAlgebraic34/45.0000);
  const double_v intermediateAlgebraic51 = intermediateAlgebraic36+intermediateAlgebraic41+intermediateAlgebraic43+intermediateAlgebraic46+intermediateAlgebraic50+- constant54;
  const double_v intermediateRate0 = - ((intermediateAlgebraic51+intermediateAlgebraic0)/constant0);
  const double_v intermediateAlgebraic32 =  intermediateState1*((intermediateState3 -  intermediateState2*exponential(( -1.00000*constant6*intermediateState1)/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*intermediateState1)/( constant35*constant36))));
  const double_v intermediateAlgebraic26 =  (( constant35*constant36)/constant6)*log(intermediateState2/intermediateState3);
  const double_v intermediateAlgebraic56 =  intermediateState2*exponential( ( - constant41*intermediateAlgebraic26)*(constant6/( constant35*constant36)));
  const double_v intermediateAlgebraic57 =  constant40*(sqr(intermediateAlgebraic56)/(constant42+sqr(intermediateAlgebraic56)));
  const double_v intermediateAlgebraic58 = 1.00000 - 1.0/(1.00000+( constant43*(1.00000+sqr(intermediateAlgebraic56)/constant42))/( pow(constant46, 2.00000)*exponential(( 2.00000*(1.00000 - constant41)*intermediateState1*constant6)/( constant35*constant36))));
  const double_v intermediateAlgebraic59 =  intermediateAlgebraic57*intermediateAlgebraic58;
  const double_v intermediateAlgebraic60 =  constant50*intermediateAlgebraic59*(intermediateAlgebraic32/50.0000);
  const double_v intermediateAlgebraic61 =  ( constant38*pow4(intermediateState13))*intermediateState14;
  const double_v intermediateAlgebraic62 =  constant51*intermediateAlgebraic61*(intermediateAlgebraic32/50.0000);
  const double_v intermediateAlgebraic66 =  (1.00000/7.00000)*(exponential(intermediateState6/67.3000) - 1.00000);
  const double_v intermediateAlgebraic67 = 1.0/(1.00000+ 0.120000*exponential( -0.100000*intermediateState1*(constant6/( constant35*constant36)))+ 0.0400000*intermediateAlgebraic66*exponential(- ( intermediateState1*(constant6/( constant35*constant36)))));
  const double_v intermediateAlgebraic68 =  constant6*(constant47/( sqr(1.00000+constant44/intermediateState2)*pow3(1.00000+constant45/intermediateState5)));
  const double_v intermediateAlgebraic69 =  constant53*intermediateAlgebraic68*intermediateAlgebraic67;
  const double_v intermediateRate3 =  - constant9*((intermediateAlgebraic60+intermediateAlgebraic62+constant12+ - 2.00000*intermediateAlgebraic69)/( (1000.00/1.00000)*constant6*constant3)) - (intermediateAlgebraic41+intermediateAlgebraic43+constant12+ -2.00000*intermediateAlgebraic50)/( (1000.00/1.00000)*constant6*constant4);
  const double_v intermediateRate2 = (intermediateAlgebraic60+intermediateAlgebraic62+constant12+ - 2.00000*intermediateAlgebraic69)/( (1000.00/1.00000)*constant6*constant3) - (intermediateState2 - intermediateState4)/constant7;
  const double_v intermediateAlgebraic64 =  ( ( constant39*pow3(intermediateState15))*intermediateState16)*intermediateState17;
  const double_v intermediateAlgebraic63 =  intermediateState1*((intermediateState5 -  intermediateState6*exponential(( -1.00000*constant6*intermediateState1)/( constant35*constant36)))/(1.00000 - exponential(( -1.00000*constant6*intermediateState1)/( constant35*constant36))));
  const double_v intermediateAlgebraic65 =  constant52*intermediateAlgebraic64*(intermediateAlgebraic63/75.0000);
  const double_v intermediateRate5 =  - constant9*((intermediateAlgebraic65+constant13+ 3.00000*intermediateAlgebraic69)/( (1000.00/1.00000)*constant6*constant3)) - (intermediateAlgebraic46+constant13+ 3.00000*intermediateAlgebraic50)/( (1000.00/1.00000)*constant6*constant4);
  const double_v intermediateRate6 = (intermediateAlgebraic65+constant13+ 3.00000*intermediateAlgebraic69)/( (1000.00/1.00000)*constant6*constant3) - (intermediateState6 - intermediateState7)/constant8;
  const double_v intermediateAlgebraic29 = 156.500/(5.00000+exponential(( - constant6*intermediateAlgebraic26)/( constant35*constant36)));
  const double_v intermediateAlgebraic30 = 156.500 -  5.00000*intermediateAlgebraic29;
  const double_v intermediateAlgebraic53 =  intermediateState1*((intermediateAlgebraic29 -  intermediateAlgebraic30*exponential(( constant6*intermediateState1)/( constant35*constant36)))/(1.00000 - exponential(( constant6*intermediateState1)/( constant35*constant36))));
  const double_v intermediateAlgebraic52 = 1.00000/(1.00000+exponential((intermediateState1 - constant23)/constant26));
  const double_v intermediateAlgebraic54 =  constant37*pow4(intermediateAlgebraic52);
  const double_v intermediateAlgebraic55 =  constant49*intermediateAlgebraic54*(intermediateAlgebraic53/45.0000);
  const double_v intermediateAlgebraic70 = intermediateAlgebraic55+intermediateAlgebraic60+intermediateAlgebraic62+intermediateAlgebraic65+intermediateAlgebraic69;
  const double_v intermediateRate1 = - ((intermediateAlgebraic70 - intermediateAlgebraic0/constant1)/constant0);


  // final step
  // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]
  states[0] += 0.5*timeStepWidth*(rate0 + intermediateRate0);
  states[1] += 0.5*timeStepWidth*(rate1 + intermediateRate1);
  states[2] += 0.5*timeStepWidth*(rate2 + intermediateRate2);
  states[3] += 0.5*timeStepWidth*(rate3 + intermediateRate3);
  states[4] += 0.5*timeStepWidth*(rate4 + intermediateRate4);
  states[5] += 0.5*timeStepWidth*(rate5 + intermediateRate5);
  states[6] += 0.5*timeStepWidth*(rate6 + intermediateRate6);
  states[7] += 0.5*timeStepWidth*(rate7 + intermediateRate7);
  states[8] += 0.5*timeStepWidth*(rate8 + intermediateRate8);
  states[9] += 0.5*timeStepWidth*(rate9 + intermediateRate9);
  states[10] += 0.5*timeStepWidth*(rate10 + intermediateRate10);
  states[11] += 0.5*timeStepWidth*(rate11 + intermediateRate11);
  states[12] += 0.5*timeStepWidth*(rate12 + intermediateRate12);
  states[13] += 0.5*timeStepWidth*(rate13 + intermediateRate13);
  states[14] += 0.5*timeStepWidth*(rate14 + intermediateRate14);
  states[15] += 0.5*timeStepWidth*(rate15 + intermediateRate15);
  states[16] += 0.5*timeStepWidth*(rate16 + intermediateRate16);
  states[17] += 0.5*timeStepWidth*(rate17 + intermediateRate17);
  states[18] += 0.5*timeStepWidth*(rate18 + intermediateRate18);
  states[19] += 0.5*timeStepWidth*(rate19 + intermediateRate19);
  states[20] += 0.5*timeStepWidth*(rate20 + intermediateRate20);
  states[21] += 0.5*timeStepWidth*(rate21 + intermediateRate21);
  states[22] += 0.5*timeStepWidth*(rate22 + intermediateRate22);
  states[23] += 0.5*timeStepWidth*(rate23 + intermediateRate23);
  states[24] += 0.5*timeStepWidth*(rate24 + intermediateRate24);
  states[25] += 0.5*timeStepWidth*(rate25 + intermediateRate25);
  states[26] += 0.5*timeStepWidth*(rate26 + intermediateRate26);
  states[27] += 0.5*timeStepWidth*(rate27 + intermediateRate27);
  states[28] += 0.5*timeStepWidth*(rate28 + intermediateRate28);
  states[29] += 0.5*timeStepWidth*(rate29 + intermediateRate29);
  states[30] += 0.5*timeStepWidth*(rate30 + intermediateRate30);
  states[31] += 0.5*timeStepWidth*(rate31 + intermediateRate31);
  states[32] += 0.5*timeStepWidth*(rate32 + intermediateRate32);
  states[33] += 0.5*timeStepWidth*(rate33 + intermediateRate33);
  states[34] += 0.5*timeStepWidth*(rate34 + intermediateRate34);
  states[35] += 0.5*timeStepWidth*(rate35 + intermediateRate35);
  states[36] += 0.5*timeStepWidth*(rate36 + intermediateRate36);
  states[37] += 0.5*timeStepWidth*(rate37 + intermediateRate37);
  states[38] += 0.5*timeStepWidth*(rate38 + intermediateRate38);
  states[39] += 0.5*timeStepWidth*(rate39 + intermediateRate39);
  states[40] += 0.5*timeStepWidth*(rate40 + intermediateRate40);
  states[41] += 0.5*timeStepWidth*(rate41 + intermediateRate41);
  states[42] += 0.5*timeStepWidth*(rate42 + intermediateRate42);
  states[43] += 0.5*timeStepWidth*(rate43 + intermediateRate43);
  states[44] += 0.5*timeStepWidth*(rate44 + intermediateRate44);
  states[45] += 0.5*timeStepWidth*(rate45 + intermediateRate45);
  states[46] += 0.5*timeStepWidth*(rate46 + intermediateRate46);
  states[47] += 0.5*timeStepWidth*(rate47 + intermediateRate47);
  states[48] += 0.5*timeStepWidth*(rate48 + intermediateRate48);
  states[49] += 0.5*timeStepWidth*(rate49 + intermediateRate49);
  states[50] += 0.5*timeStepWidth*(rate50 + intermediateRate50);
  states[51] += 0.5*timeStepWidth*(rate51 + intermediateRate51);
  states[52] += 0.5*timeStepWidth*(rate52 + intermediateRate52);
  states[53] += 0.5*timeStepWidth*(rate53 + intermediateRate53);
  states[54] += 0.5*timeStepWidth*(rate54 + intermediateRate54);
  states[55] += 0.5*timeStepWidth*(rate55 + intermediateRate55);
  states[56] += 0.5*timeStepWidth*(rate56 + intermediateRate56);

  if (stimulate)
  {
    for (int i = 0; i < std::min(3,(int)Vc::double_v::Size); i++)
    {
      states[0][i] = 20.0;
    }
  }
  //LOG(INFO) << "min: " << minX << ", max: " << maxX;
}

//! set the initial values for all states
void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<57, 71, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
initializeStates(Vc::double_v states[])
{
  states[0] = -79.974;
  states[1] = -80.2;
  states[2] = 5.9;
  states[3] = 150.9;
  states[4] = 5.9;
  states[5] = 12.7;
  states[6] = 132.0;
  states[7] = 133.0;
  states[8] = 0.009466;
  states[9] = 0.9952;
  states[10] = 0.0358;
  states[11] = 0.4981;
  states[12] = 0.581;
  states[13] = 0.009466;
  states[14] = 0.9952;
  states[15] = 0.0358;
  states[16] = 0.4981;
  states[17] = 0.581;
  states[18] = 0.0;
  states[19] = 0.0;
  states[20] = 0.0;
  states[21] = 0.0;
  states[22] = 0.0;
  states[23] = 1.0;
  states[24] = 0.0;
  states[25] = 0.0;
  states[26] = 0.0;
  states[27] = 0.0;
  states[28] = 0.0;
  states[29] = 0.1;
  states[30] = 1500.0;
  states[31] = 0.1;
  states[32] = 1500.0;
  states[33] = 25;
  states[34] = 615.000000;
  states[35] = 615.000000;
  states[36] = 811.000000;
  states[37] = 811.000000;
  states[38] = 16900.0;
  states[39] = 16900.0;
  states[40] = 0.4;
  states[41] = 0.4;
  states[42] = 7200.0;
  states[43] = 7200.0;
  states[44] = 799.6;
  states[45] = 799.6;
  states[46] = 1000.0;
  states[47] = 1000.0;
  states[48] = 3.0;
  states[49] = 0.8;
  states[50] = 1.2;
  states[51] = 3.0;
  states[52] = 0.3;
  states[53] = 0.23;
  states[54] = 0.23;
  states[55] = 0.23;
  states[56] = 0.23;
}