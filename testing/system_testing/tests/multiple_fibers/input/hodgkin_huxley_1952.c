/*
   There are a total of 9 entries in the algebraic variable array.
   There are a total of 4 entries in each of the rate and state variable arrays.
   There are a total of 9 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is E_R in component membrane (millivolt).
 * CONSTANTS[1] is Cm in component membrane (microF_per_cm2).
 * ALGEBRAIC[0] is i_Na in component sodium_channel (microA_per_cm2).
 * ALGEBRAIC[4] is i_K in component potassium_channel (microA_per_cm2).
 * ALGEBRAIC[8] is i_L in component leakage_current (microA_per_cm2).
 * CONSTANTS[2] is i_Stim in component membrane (microA_per_cm2).
 * CONSTANTS[3] is g_Na in component sodium_channel (milliS_per_cm2).
 * CONSTANTS[6] is E_Na in component sodium_channel (millivolt).
 * STATES[1] is m in component sodium_channel_m_gate (dimensionless).
 * STATES[2] is h in component sodium_channel_h_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component sodium_channel_m_gate (per_millisecond).
 * ALGEBRAIC[5] is beta_m in component sodium_channel_m_gate (per_millisecond).
 * ALGEBRAIC[2] is alpha_h in component sodium_channel_h_gate (per_millisecond).
 * ALGEBRAIC[6] is beta_h in component sodium_channel_h_gate (per_millisecond).
 * CONSTANTS[4] is g_K in component potassium_channel (milliS_per_cm2).
 * CONSTANTS[7] is E_K in component potassium_channel (millivolt).
 * STATES[3] is n in component potassium_channel_n_gate (dimensionless).
 * ALGEBRAIC[3] is alpha_n in component potassium_channel_n_gate (per_millisecond).
 * ALGEBRAIC[7] is beta_n in component potassium_channel_n_gate (per_millisecond).
 * CONSTANTS[5] is g_L in component leakage_current (milliS_per_cm2).
 * CONSTANTS[8] is E_L in component leakage_current (millivolt).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[1] is d/dt m in component sodium_channel_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_channel_h_gate (dimensionless).
 * RATES[3] is d/dt n in component potassium_channel_n_gate (dimensionless).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = -75;
CONSTANTS[0] = -75;
CONSTANTS[1] = 1;
CONSTANTS[2] = 0;
CONSTANTS[3] = 120;
STATES[1] = 0.05;
STATES[2] = 0.6;
CONSTANTS[4] = 36;
STATES[3] = 0.325;
CONSTANTS[5] = 0.3;
CONSTANTS[6] = CONSTANTS[0]+115.000;
CONSTANTS[7] = CONSTANTS[0] - 12.0000;
CONSTANTS[8] = CONSTANTS[0]+10.6130;
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[1] = ( - 0.100000*(STATES[0]+50.0000))/(exp(- (STATES[0]+50.0000)/10.0000) - 1.00000);
ALGEBRAIC[5] =  4.00000*exp(- (STATES[0]+75.0000)/18.0000);
RATES[1] =  ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[5]*STATES[1];
ALGEBRAIC[2] =  0.0700000*exp(- (STATES[0]+75.0000)/20.0000);
ALGEBRAIC[6] = 1.00000/(exp(- (STATES[0]+45.0000)/10.0000)+1.00000);
RATES[2] =  ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[6]*STATES[2];
ALGEBRAIC[3] = ( - 0.0100000*(STATES[0]+65.0000))/(exp(- (STATES[0]+65.0000)/10.0000) - 1.00000);
ALGEBRAIC[7] =  0.125000*exp((STATES[0]+75.0000)/80.0000);
RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[7]*STATES[3];
ALGEBRAIC[0] =  CONSTANTS[3]*pow(STATES[1], 3.00000)*STATES[2]*(STATES[0] - CONSTANTS[6]);
ALGEBRAIC[4] =  CONSTANTS[4]*pow(STATES[3], 4.00000)*(STATES[0] - CONSTANTS[7]);
ALGEBRAIC[8] =  CONSTANTS[5]*(STATES[0] - CONSTANTS[8]);
RATES[0] = - (- CONSTANTS[2]+ALGEBRAIC[0]+ALGEBRAIC[4]+ALGEBRAIC[8])/CONSTANTS[1];
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[1] = ( - 0.100000*(STATES[0]+50.0000))/(exp(- (STATES[0]+50.0000)/10.0000) - 1.00000);
ALGEBRAIC[5] =  4.00000*exp(- (STATES[0]+75.0000)/18.0000);
ALGEBRAIC[2] =  0.0700000*exp(- (STATES[0]+75.0000)/20.0000);
ALGEBRAIC[6] = 1.00000/(exp(- (STATES[0]+45.0000)/10.0000)+1.00000);
ALGEBRAIC[3] = ( - 0.0100000*(STATES[0]+65.0000))/(exp(- (STATES[0]+65.0000)/10.0000) - 1.00000);
ALGEBRAIC[7] =  0.125000*exp((STATES[0]+75.0000)/80.0000);
ALGEBRAIC[0] =  CONSTANTS[3]*pow(STATES[1], 3.00000)*STATES[2]*(STATES[0] - CONSTANTS[6]);
ALGEBRAIC[4] =  CONSTANTS[4]*pow(STATES[3], 4.00000)*(STATES[0] - CONSTANTS[7]);
ALGEBRAIC[8] =  CONSTANTS[5]*(STATES[0] - CONSTANTS[8]);
}
