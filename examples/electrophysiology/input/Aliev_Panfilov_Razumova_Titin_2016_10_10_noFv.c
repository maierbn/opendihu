/*
   There are a total of 8 entries in the algebraic variable array.
   There are a total of 9 entries in each of the rate and state variable arrays.
   There are a total of 38 entries in the constant variable array.
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
 * CONSTANTS[6] is stim in component Aliev_Panfilov (dimensionless).
 * STATES[2] is D in component Razumova (dimensionless).
 * STATES[3] is A_1 in component Razumova (dimensionless).
 * STATES[4] is A_2 in component Razumova (dimensionless).
 * STATES[5] is x_1 in component Razumova (micrometer).
 * STATES[6] is x_2 in component Razumova (micrometer).
 * CONSTANTS[7] is x_0 in component Razumova (micrometer).
 * CONSTANTS[8] is R_T in component Razumova (dimensionless).
 * ALGEBRAIC[7] is R_off in component Razumova (dimensionless).
 * CONSTANTS[9] is A_2_0 in component Razumova (dimensionless).
 * CONSTANTS[10] is r_0 in component Razumova (dimensionless).
 * CONSTANTS[11] is l_hs in component Razumova (micrometer).
 * CONSTANTS[12] is rel_velo in component Razumova (micrometer_per_millisecond).
 * CONSTANTS[30] is velo_scaled in component Razumova (micrometer_per_millisecond).
 * CONSTANTS[13] is velo_max in component Razumova (dimensionless).
 * ALGEBRAIC[0] is k_on in component Razumova (per_millisecond).
 * ALGEBRAIC[6] is k_off in component Razumova (per_millisecond).
 * CONSTANTS[31] is f in component Razumova (per_millisecond).
 * CONSTANTS[32] is f_prime in component Razumova (per_millisecond).
 * CONSTANTS[33] is h in component Razumova (per_millisecond).
 * CONSTANTS[34] is h_prime in component Razumova (per_millisecond).
 * CONSTANTS[35] is g in component Razumova (per_millisecond).
 * CONSTANTS[36] is g_prime in component Razumova (per_millisecond).
 * CONSTANTS[14] is k_0_on in component Razumova (per_millisecond).
 * CONSTANTS[15] is k_0_off in component Razumova (per_millisecond).
 * CONSTANTS[16] is k_Ca_on in component Razumova (per_millisecond).
 * CONSTANTS[17] is k_Ca_off in component Razumova (per_millisecond).
 * CONSTANTS[18] is g_0 in component Razumova (per_millisecond).
 * CONSTANTS[19] is f_0 in component Razumova (per_millisecond).
 * CONSTANTS[20] is h_0 in component Razumova (per_millisecond).
 * CONSTANTS[21] is f_prime0 in component Razumova (per_millisecond).
 * CONSTANTS[22] is h_prime0 in component Razumova (per_millisecond).
 * CONSTANTS[23] is g_prime0 in component Razumova (per_millisecond).
 * ALGEBRAIC[1] is sigma in component Razumova (dimensionless).
 * ALGEBRAIC[2] is lambda_A1 in component Razumova (dimensionless).
 * ALGEBRAIC[3] is lambda_A2 in component Razumova (dimensionless).
 * CONSTANTS[24] is Ca_50 in component Razumova (dimensionless).
 * CONSTANTS[25] is nu in component Razumova (dimensionless).
 * CONSTANTS[26] is E_ATP in component Razumova (joule).
 * CONSTANTS[27] is kappa in component Razumova (joule_per_kelvin).
 * CONSTANTS[28] is T in component Razumova (kelvin).
 * ALGEBRAIC[4] is ActiveStress in component Razumova (dimensionless).
 * ALGEBRAIC[5] is Activation in component Razumova (dimensionless).
 * CONSTANTS[37] is f_l in component Razumova (dimensionless).
 * CONSTANTS[29] is A_2_max in component Razumova (dimensionless).
 * STATES[7] is Dummy_1 in component Razumova (dimensionless).
 * STATES[8] is Dummy_2 in component Razumova (dimensionless).
 * RATES[0] is d/dt V_m in component Aliev_Panfilov (dimensionless).
 * RATES[1] is d/dt r in component Aliev_Panfilov (dimensionless).
 * RATES[2] is d/dt D in component Razumova (dimensionless).
 * RATES[3] is d/dt A_1 in component Razumova (dimensionless).
 * RATES[4] is d/dt A_2 in component Razumova (dimensionless).
 * RATES[5] is d/dt x_1 in component Razumova (micrometer).
 * RATES[6] is d/dt x_2 in component Razumova (micrometer).
 * RATES[7] is d/dt Dummy_1 in component Razumova (dimensionless).
 * RATES[8] is d/dt Dummy_2 in component Razumova (dimensionless).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
STATES[0] = 0;
STATES[1] = 0;
CONSTANTS[0] = 0;
CONSTANTS[1] = 128;
CONSTANTS[2] = 0.15;
CONSTANTS[3] = 0.002;
CONSTANTS[4] = 0.2;
CONSTANTS[5] = 0.3;
CONSTANTS[6] = 1;
STATES[2] = 3.8e-14;
STATES[3] = 1e-14;
STATES[4] = 3.4e-13;
STATES[5] = 1e-16;
STATES[6] = 8e-3;
CONSTANTS[7] = 8e-3;
CONSTANTS[8] = 1;
CONSTANTS[9] = 3.4e-13;
CONSTANTS[10] = 3.0722;
CONSTANTS[11] = 1;
CONSTANTS[12] = 0;
CONSTANTS[13] = 7.815e-5;
CONSTANTS[14] = 0;
CONSTANTS[15] = 100e-3;
CONSTANTS[16] = 120e-3;
CONSTANTS[17] = 50e-3;
CONSTANTS[18] = 4e-3;
CONSTANTS[19] = 50e-3;
CONSTANTS[20] = 8e-3;
CONSTANTS[21] = 400e-3;
CONSTANTS[22] = 6e-3;
CONSTANTS[23] = 3.5400e-13;
CONSTANTS[24] = 1;
CONSTANTS[25] = 3.2;
CONSTANTS[26] = 9.1362e-20;
CONSTANTS[27] = 1.38e-23;
CONSTANTS[28] = 310;
CONSTANTS[29] = 0.0444557705319849;
STATES[7] = 0;
STATES[8] = 0;
CONSTANTS[30] =  CONSTANTS[12]*CONSTANTS[13]*0.00000;
CONSTANTS[31] = CONSTANTS[19];
CONSTANTS[32] = CONSTANTS[21];
CONSTANTS[33] = CONSTANTS[20];
CONSTANTS[34] = CONSTANTS[22];
CONSTANTS[35] = CONSTANTS[18];
CONSTANTS[36] = 0.00000;
CONSTANTS[37] = (CONSTANTS[11]<0.635000 ? 0.00000 : CONSTANTS[11]<0.835000 ?  4.20000*(CONSTANTS[11] - 0.635000) : CONSTANTS[11]<1.00000 ? 0.840000+ 0.969700*(CONSTANTS[11] - 0.835000) : CONSTANTS[11]<1.12500 ? 1.00000 : CONSTANTS[11]<1.82500 ? 1.00000 -  1.42860*(CONSTANTS[11] - 1.12500) : 0.00000);
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[0] = ( - CONSTANTS[1]*STATES[0]*(STATES[0] - CONSTANTS[2])*(STATES[0] - 1.00000) -  STATES[0]*STATES[1])+CONSTANTS[0];
RATES[1] =  (CONSTANTS[3]+( CONSTANTS[4]*STATES[1])/(CONSTANTS[5]+STATES[0]))*(- STATES[1] -  CONSTANTS[1]*STATES[0]*((STATES[0] - CONSTANTS[2]) - 1.00000));
RATES[3] = ( CONSTANTS[31]*STATES[2]+ CONSTANTS[34]*STATES[4]) -  (CONSTANTS[32]+CONSTANTS[33])*STATES[3];
RATES[4] = ( CONSTANTS[33]*STATES[3] -  (CONSTANTS[34]+CONSTANTS[35])*STATES[4])+ CONSTANTS[36]*STATES[2];
RATES[5] = ( (( - CONSTANTS[31]*STATES[2])/STATES[3])*STATES[5] -  (( CONSTANTS[34]*STATES[4])/STATES[3])*STATES[5])+CONSTANTS[30];
RATES[6] =  (( - CONSTANTS[33]*STATES[3])/STATES[4])*(STATES[6] - CONSTANTS[7])+CONSTANTS[30];
ALGEBRAIC[4] =  (( 1.12500*(( STATES[4]*STATES[6]+ STATES[3]*STATES[5]) -  CONSTANTS[9]*CONSTANTS[7]))/( CONSTANTS[29]*CONSTANTS[7]))*CONSTANTS[37];
RATES[7] = ALGEBRAIC[4];
ALGEBRAIC[5] = (STATES[4] - CONSTANTS[9])/CONSTANTS[29];
RATES[8] = ALGEBRAIC[5];
ALGEBRAIC[7] = ((CONSTANTS[8] - STATES[3]) - STATES[4]) - STATES[2];
ALGEBRAIC[0] = CONSTANTS[14]+( (CONSTANTS[16] - CONSTANTS[14])*STATES[1]*CONSTANTS[10])/( STATES[1]*CONSTANTS[10]+CONSTANTS[24]);
ALGEBRAIC[6] = CONSTANTS[15]+( (CONSTANTS[17] - CONSTANTS[15])*STATES[1]*CONSTANTS[10])/( STATES[1]*CONSTANTS[10]+CONSTANTS[24]);
RATES[2] = ( ALGEBRAIC[0]*ALGEBRAIC[7]+ CONSTANTS[32]*STATES[3]+ CONSTANTS[35]*STATES[4]) -  (ALGEBRAIC[6]+CONSTANTS[31]+CONSTANTS[36])*STATES[2];
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[4] =  (( 1.12500*(( STATES[4]*STATES[6]+ STATES[3]*STATES[5]) -  CONSTANTS[9]*CONSTANTS[7]))/( CONSTANTS[29]*CONSTANTS[7]))*CONSTANTS[37];
ALGEBRAIC[5] = (STATES[4] - CONSTANTS[9])/CONSTANTS[29];
ALGEBRAIC[7] = ((CONSTANTS[8] - STATES[3]) - STATES[4]) - STATES[2];
ALGEBRAIC[0] = CONSTANTS[14]+( (CONSTANTS[16] - CONSTANTS[14])*STATES[1]*CONSTANTS[10])/( STATES[1]*CONSTANTS[10]+CONSTANTS[24]);
ALGEBRAIC[6] = CONSTANTS[15]+( (CONSTANTS[17] - CONSTANTS[15])*STATES[1]*CONSTANTS[10])/( STATES[1]*CONSTANTS[10]+CONSTANTS[24]);
ALGEBRAIC[1] = (STATES[6]>CONSTANTS[7] ? 1.00000 : STATES[6]<CONSTANTS[7] ? 8.00000 : 0.00000);
ALGEBRAIC[2] = STATES[3]/CONSTANTS[8];
ALGEBRAIC[3] = STATES[4]/CONSTANTS[8];
}
