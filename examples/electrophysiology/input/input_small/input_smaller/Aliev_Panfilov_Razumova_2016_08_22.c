/*
   There are a total of 12 entries in the algebraic variable array.
   There are a total of 7 entries in each of the rate and state variable arrays.
   There are a total of 30 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V_m in component Aliev_Panfilov (dimensionless).
 * STATES[1] is r in component Aliev_Panfilov (dimensionless).
 * CONSTANTS[0] is I_HH in component Aliev_Panfilov (dimensionless).
 * CONSTANTS[1] is k in component Aliev_Panfilov (dimensionless).
 * CONSTANTS[2] is a in component Aliev_Panfilov (dimensionless).
 * CONSTANTS[3] is e_0 in component Aliev_Panfilov (dimensionless).
 * CONSTANTS[4] is m_1 in component Aliev_Panfilov (dimensionless).
 * CONSTANTS[5] is m_2 in component Aliev_Panfilov (dimensionless).
 * STATES[2] is D in component Razumova (dimensionless).
 * STATES[3] is A_1 in component Razumova (dimensionless).
 * STATES[4] is A_2 in component Razumova (dimensionless).
 * STATES[5] is x_1 in component Razumova (micrometer).
 * STATES[6] is x_2 in component Razumova (micrometer).
 * CONSTANTS[6] is x_0 in component Razumova (micrometer).
 * CONSTANTS[7] is R_T in component Razumova (dimensionless).
 * ALGEBRAIC[11] is R_off in component Razumova (dimensionless).
 * CONSTANTS[8] is l_hs in component Razumova (micrometer).
 * CONSTANTS[9] is velo in component Razumova (micrometer_per_millisecond).
 * CONSTANTS[29] is velo_scaled in component Razumova (micrometer_per_millisecond).
 * CONSTANTS[10] is velo_max in component Razumova (dimensionless).
 * ALGEBRAIC[3] is k_on in component Razumova (per_millisecond).
 * ALGEBRAIC[4] is k_off in component Razumova (per_millisecond).
 * ALGEBRAIC[5] is f in component Razumova (per_millisecond).
 * ALGEBRAIC[6] is f_prime in component Razumova (per_millisecond).
 * ALGEBRAIC[7] is h in component Razumova (per_millisecond).
 * ALGEBRAIC[8] is h_prime in component Razumova (per_millisecond).
 * ALGEBRAIC[9] is g in component Razumova (per_millisecond).
 * ALGEBRAIC[10] is g_prime in component Razumova (per_millisecond).
 * CONSTANTS[11] is k_0_on in component Razumova (per_millisecond).
 * CONSTANTS[12] is k_0_off in component Razumova (per_millisecond).
 * CONSTANTS[13] is k_Ca_on in component Razumova (per_millisecond).
 * CONSTANTS[14] is k_Ca_off in component Razumova (per_millisecond).
 * CONSTANTS[15] is g_0 in component Razumova (per_millisecond).
 * CONSTANTS[16] is f_0 in component Razumova (per_millisecond).
 * CONSTANTS[17] is h_0 in component Razumova (per_millisecond).
 * CONSTANTS[18] is f_prime0 in component Razumova (per_millisecond).
 * CONSTANTS[19] is h_prime0 in component Razumova (per_millisecond).
 * CONSTANTS[20] is g_prime0 in component Razumova (per_millisecond).
 * ALGEBRAIC[0] is sigma in component Razumova (dimensionless).
 * ALGEBRAIC[1] is lambda_A1 in component Razumova (dimensionless).
 * ALGEBRAIC[2] is lambda_A2 in component Razumova (dimensionless).
 * CONSTANTS[21] is Ca_50 in component Razumova (dimensionless).
 * CONSTANTS[22] is nu in component Razumova (dimensionless).
 * CONSTANTS[23] is E_ATP in component Razumova (joule).
 * CONSTANTS[24] is kappa in component Razumova (joule_per_kelvin).
 * CONSTANTS[25] is T in component Razumova (kelvin).
 * CONSTANTS[26] is k_on_rel in component Razumova (dimensionless).
 * CONSTANTS[27] is k_off_rel in component Razumova (dimensionless).
 * CONSTANTS[28] is alpha in component Razumova (dimensionless).
 * RATES[0] is d/dt V_m in component Aliev_Panfilov (dimensionless).
 * RATES[1] is d/dt r in component Aliev_Panfilov (dimensionless).
 * RATES[2] is d/dt D in component Razumova (dimensionless).
 * RATES[3] is d/dt A_1 in component Razumova (dimensionless).
 * RATES[4] is d/dt A_2 in component Razumova (dimensionless).
 * RATES[5] is d/dt x_1 in component Razumova (micrometer).
 * RATES[6] is d/dt x_2 in component Razumova (micrometer).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = 0;
STATES[1] = 20;
CONSTANTS[0] = 0.0;
CONSTANTS[1] = 128;
CONSTANTS[2] = 0.15;
CONSTANTS[3] = 0.002;
CONSTANTS[4] = 0.2;
CONSTANTS[5] = 0.3;
STATES[2] = 0.0001;
STATES[3] = 0.0001;
STATES[4] = 0.0001;
STATES[5] = 1e-16;
STATES[6] = 8e-3;
CONSTANTS[6] = 8e-3;
CONSTANTS[7] = 1;
CONSTANTS[8] = 1;
CONSTANTS[9] = 0;
CONSTANTS[10] = 1e-4;
CONSTANTS[11] = 0;
CONSTANTS[12] = 100e-3;
CONSTANTS[13] = 120e-3;
CONSTANTS[14] = 50e-3;
CONSTANTS[15] = 4e-3;
CONSTANTS[16] = 50e-3;
CONSTANTS[17] = 8e-3;
CONSTANTS[18] = 400e-3;
CONSTANTS[19] = 6e-3;
CONSTANTS[20] = 3.5400e-13;
CONSTANTS[21] = 26;
CONSTANTS[22] = 3.2;
CONSTANTS[23] = 9.1362e-20;
CONSTANTS[24] = 1.38e-23;
CONSTANTS[25] = 310;
CONSTANTS[26] = 0.925;
CONSTANTS[27] = 0.143;
CONSTANTS[28] = 1.35;
CONSTANTS[29] =  CONSTANTS[9]*CONSTANTS[10];
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[0] = ( - CONSTANTS[1]*STATES[0]*(STATES[0] - CONSTANTS[2])*(STATES[0] - 1.00000) -  STATES[0]*STATES[1])+CONSTANTS[0];
RATES[1] =  (CONSTANTS[3]+( CONSTANTS[4]*STATES[1])/(CONSTANTS[5]+STATES[0]))*(- STATES[1] -  CONSTANTS[1]*STATES[0]*((STATES[0] - CONSTANTS[2]) - 1.00000));
ALGEBRAIC[0] = (STATES[6]>CONSTANTS[6] ? 1.00000 : STATES[6]<CONSTANTS[6] ? 8.00000 : 0.00000);
ALGEBRAIC[7] =  CONSTANTS[28]*CONSTANTS[17]*exp( ALGEBRAIC[0]*pow(STATES[5], 2.00000));
RATES[6] =  (( - ALGEBRAIC[7]*STATES[3])/STATES[4])*(STATES[6] - CONSTANTS[6])+CONSTANTS[29];
ALGEBRAIC[1] = STATES[3]/CONSTANTS[7];
ALGEBRAIC[2] = STATES[4]/CONSTANTS[7];
ALGEBRAIC[5] =  CONSTANTS[28]*CONSTANTS[16]*pow(1.00000+ ALGEBRAIC[1]*(exp( (STATES[5]/CONSTANTS[6])*(CONSTANTS[22] - 1.00000)) - 1.00000)+ ALGEBRAIC[2]*(exp( (STATES[6]/CONSTANTS[6])*(CONSTANTS[22] - 1.00000)) - 1.00000), 2.00000);
ALGEBRAIC[6] =  CONSTANTS[28]*CONSTANTS[18]*exp( ALGEBRAIC[0]*pow(STATES[5], 2.00000));
ALGEBRAIC[8] =  CONSTANTS[28]*CONSTANTS[19]*exp( ALGEBRAIC[0]*pow(STATES[6], 2.00000) -  ALGEBRAIC[0]*pow(STATES[5], 2.00000));
RATES[3] = ( ALGEBRAIC[5]*STATES[2]+ ALGEBRAIC[8]*STATES[4]) -  (ALGEBRAIC[6]+ALGEBRAIC[7])*STATES[3];
RATES[5] = ( (( - ALGEBRAIC[5]*STATES[2])/STATES[3])*STATES[5] -  (( ALGEBRAIC[8]*STATES[4])/STATES[3])*STATES[5])+CONSTANTS[29];
ALGEBRAIC[9] =  CONSTANTS[28]*CONSTANTS[15]*exp( ALGEBRAIC[0]*pow(STATES[6] - CONSTANTS[6], 2.00000));
ALGEBRAIC[10] =  (( ALGEBRAIC[5]*ALGEBRAIC[7]*ALGEBRAIC[9])/( ALGEBRAIC[6]*ALGEBRAIC[8]))*exp(- CONSTANTS[23]/( CONSTANTS[24]*CONSTANTS[25]));
RATES[4] = ( ALGEBRAIC[7]*STATES[3] -  (ALGEBRAIC[8]+ALGEBRAIC[9])*STATES[4])+ ALGEBRAIC[10]*STATES[2];
ALGEBRAIC[11] = ((CONSTANTS[7] - STATES[3]) - STATES[4]) - STATES[2];
ALGEBRAIC[3] =  CONSTANTS[26]*(CONSTANTS[11]+( (CONSTANTS[13] - CONSTANTS[11])*STATES[1])/(STATES[1]+CONSTANTS[21]));
ALGEBRAIC[4] =  CONSTANTS[27]*(CONSTANTS[12]+( (CONSTANTS[14] - CONSTANTS[12])*STATES[1])/(STATES[1]+CONSTANTS[21]));
RATES[2] = ( ALGEBRAIC[3]*ALGEBRAIC[11]+ ALGEBRAIC[6]*STATES[3]+ ALGEBRAIC[9]*STATES[4]) -  (ALGEBRAIC[4]+ALGEBRAIC[5]+ALGEBRAIC[10])*STATES[2];
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[0] = (STATES[6]>CONSTANTS[6] ? 1.00000 : STATES[6]<CONSTANTS[6] ? 8.00000 : 0.00000);
ALGEBRAIC[7] =  CONSTANTS[28]*CONSTANTS[17]*exp( ALGEBRAIC[0]*pow(STATES[5], 2.00000));
ALGEBRAIC[1] = STATES[3]/CONSTANTS[7];
ALGEBRAIC[2] = STATES[4]/CONSTANTS[7];
ALGEBRAIC[5] =  CONSTANTS[28]*CONSTANTS[16]*pow(1.00000+ ALGEBRAIC[1]*(exp( (STATES[5]/CONSTANTS[6])*(CONSTANTS[22] - 1.00000)) - 1.00000)+ ALGEBRAIC[2]*(exp( (STATES[6]/CONSTANTS[6])*(CONSTANTS[22] - 1.00000)) - 1.00000), 2.00000);
ALGEBRAIC[6] =  CONSTANTS[28]*CONSTANTS[18]*exp( ALGEBRAIC[0]*pow(STATES[5], 2.00000));
ALGEBRAIC[8] =  CONSTANTS[28]*CONSTANTS[19]*exp( ALGEBRAIC[0]*pow(STATES[6], 2.00000) -  ALGEBRAIC[0]*pow(STATES[5], 2.00000));
ALGEBRAIC[9] =  CONSTANTS[28]*CONSTANTS[15]*exp( ALGEBRAIC[0]*pow(STATES[6] - CONSTANTS[6], 2.00000));
ALGEBRAIC[10] =  (( ALGEBRAIC[5]*ALGEBRAIC[7]*ALGEBRAIC[9])/( ALGEBRAIC[6]*ALGEBRAIC[8]))*exp(- CONSTANTS[23]/( CONSTANTS[24]*CONSTANTS[25]));
ALGEBRAIC[11] = ((CONSTANTS[7] - STATES[3]) - STATES[4]) - STATES[2];
ALGEBRAIC[3] =  CONSTANTS[26]*(CONSTANTS[11]+( (CONSTANTS[13] - CONSTANTS[11])*STATES[1])/(STATES[1]+CONSTANTS[21]));
ALGEBRAIC[4] =  CONSTANTS[27]*(CONSTANTS[12]+( (CONSTANTS[14] - CONSTANTS[12])*STATES[1])/(STATES[1]+CONSTANTS[21]));
}
