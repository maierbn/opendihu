/*
   There are a total of 71 entries in the algebraic variable array.
   There are a total of 56 entries in each of the rate and state variable arrays.
   There are a total of 105 entries in the constant variable array.
 */
/*
 * VOI is time in component cell (millisecond).
 * CONSTANTS[0] is C_m in component wal_environment (microF_per_cm2).
 * CONSTANTS[1] is gam in component wal_environment (dimensionless).
 * CONSTANTS[2] is R_a in component wal_environment (ohm_cm2).
 * CONSTANTS[3] is tsi in component wal_environment (centi_metre).
 * CONSTANTS[4] is tsi2 in component wal_environment (centi_metre).
 * CONSTANTS[5] is tsi3 in component wal_environment (centi_metre).
 * CONSTANTS[6] is FF in component wal_environment (C_per_mol).
 * CONSTANTS[7] is tau_K in component wal_environment (millisecond).
 * CONSTANTS[8] is tau_Na in component wal_environment (millisecond).
 * CONSTANTS[9] is f_T in component wal_environment (dimensionless).
 * CONSTANTS[10] is tau_K2 in component wal_environment (millisecond).
 * CONSTANTS[11] is tau_Na2 in component wal_environment (millisecond).
 * CONSTANTS[12] is I_K_rest in component wal_environment (microA_per_cm2).
 * CONSTANTS[13] is I_Na_rest in component wal_environment (microA_per_cm2).
 * CONSTANTS[14] is alpha_h_bar in component wal_environment (per_millisecond).
 * CONSTANTS[15] is alpha_m_bar in component wal_environment (per_millisecond_per_millivolt).
 * CONSTANTS[16] is alpha_n_bar in component wal_environment (per_millisecond_per_millivolt).
 * CONSTANTS[17] is beta_h_bar in component wal_environment (per_millisecond).
 * CONSTANTS[18] is beta_m_bar in component wal_environment (per_millisecond).
 * CONSTANTS[19] is beta_n_bar in component wal_environment (per_millisecond).
 * CONSTANTS[20] is V_m in component wal_environment (millivolt).
 * CONSTANTS[21] is V_n in component wal_environment (millivolt).
 * CONSTANTS[22] is V_h in component wal_environment (millivolt).
 * CONSTANTS[23] is V_a in component wal_environment (millivolt).
 * CONSTANTS[24] is V_S_inf in component wal_environment (millivolt).
 * CONSTANTS[25] is V_h_K_inf in component wal_environment (millivolt).
 * CONSTANTS[26] is A_a in component wal_environment (millivolt).
 * CONSTANTS[27] is A_S_inf in component wal_environment (millivolt).
 * CONSTANTS[28] is A_h_K_inf in component wal_environment (millivolt).
 * CONSTANTS[29] is K_alpha_h in component wal_environment (millivolt).
 * CONSTANTS[30] is K_beta_h in component wal_environment (millivolt).
 * CONSTANTS[31] is K_alpha_m in component wal_environment (millivolt).
 * CONSTANTS[32] is K_alpha_n in component wal_environment (millivolt).
 * CONSTANTS[33] is K_beta_m in component wal_environment (millivolt).
 * CONSTANTS[34] is K_beta_n in component wal_environment (millivolt).
 * CONSTANTS[35] is RR in component wal_environment (milliJ_per_degreeK_per_mol).
 * CONSTANTS[36] is TT in component wal_environment (degreeK).
 * CONSTANTS[37] is g_Cl_bar in component wal_environment (milliS_per_cm2).
 * CONSTANTS[38] is g_K_bar in component wal_environment (milliS_per_cm2).
 * CONSTANTS[39] is g_Na_bar in component wal_environment (milliS_per_cm2).
 * CONSTANTS[40] is G_K in component wal_environment (milliS_per_cm2).
 * CONSTANTS[41] is del in component wal_environment (dimensionless).
 * CONSTANTS[42] is K_K in component wal_environment (milliM2).
 * CONSTANTS[43] is K_S in component wal_environment (milliM2).
 * CONSTANTS[44] is K_m_K in component wal_environment (milliM).
 * CONSTANTS[45] is K_m_Na in component wal_environment (milliM).
 * CONSTANTS[46] is S_i in component wal_environment (milliM).
 * CONSTANTS[47] is J_NaK_bar in component wal_environment (micro_mol_per_cm2_per_second).
 * CONSTANTS[48] is V_tau in component wal_environment (millivolt).
 * ALGEBRAIC[0] is I_T in component wal_environment (microA_per_cm2).
 * STATES[0] is vS in component wal_environment (millivolt).
 * STATES[1] is vT in component wal_environment (millivolt).
 * ALGEBRAIC[51] is I_ionic_s in component wal_environment (microA_per_cm2).
 * ALGEBRAIC[70] is I_ionic_t in component wal_environment (microA_per_cm2).
 * STATES[2] is K_t in component wal_environment (milliM).
 * STATES[3] is K_i in component wal_environment (milliM).
 * STATES[4] is K_e in component wal_environment (milliM).
 * STATES[5] is Na_i in component wal_environment (milliM).
 * STATES[6] is Na_t in component wal_environment (milliM).
 * STATES[7] is Na_e in component wal_environment (milliM).
 * ALGEBRAIC[13] is E_K in component wal_environment (millivolt).
 * ALGEBRAIC[25] is E_K_t in component wal_environment (millivolt).
 * ALGEBRAIC[26] is Cl_i in component wal_environment (milliM).
 * ALGEBRAIC[27] is Cl_o in component wal_environment (milliM).
 * ALGEBRAIC[28] is Cl_i_t in component wal_environment (milliM).
 * ALGEBRAIC[29] is Cl_o_t in component wal_environment (milliM).
 * ALGEBRAIC[30] is J_K in component wal_environment (milliV_milliM).
 * ALGEBRAIC[31] is J_K_t in component wal_environment (milliV_milliM).
 * CONSTANTS[49] is eta_Cl in component wal_environment (dimensionless).
 * CONSTANTS[50] is eta_IR in component wal_environment (dimensionless).
 * CONSTANTS[51] is eta_DR in component wal_environment (dimensionless).
 * CONSTANTS[52] is eta_Na in component wal_environment (dimensionless).
 * CONSTANTS[53] is eta_NaK in component wal_environment (dimensionless).
 * ALGEBRAIC[36] is I_Cl in component sarco_Cl_channel (microA_per_cm2).
 * ALGEBRAIC[41] is I_IR in component sarco_IR_channel (microA_per_cm2).
 * ALGEBRAIC[43] is I_DR in component sarco_DR_channel (microA_per_cm2).
 * ALGEBRAIC[46] is I_Na in component sarco_Na_channel (microA_per_cm2).
 * ALGEBRAIC[50] is I_NaK in component sarco_NaK_channel (microA_per_cm2).
 * ALGEBRAIC[55] is I_Cl_t in component t_Cl_channel (microA_per_cm2).
 * ALGEBRAIC[60] is I_IR_t in component t_IR_channel (microA_per_cm2).
 * ALGEBRAIC[62] is I_DR_t in component t_DR_channel (microA_per_cm2).
 * ALGEBRAIC[65] is I_Na_t in component t_Na_channel (microA_per_cm2).
 * ALGEBRAIC[69] is I_NaK_t in component t_NaK_channel (microA_per_cm2).
 * ALGEBRAIC[32] is I_HH in component wal_environment (microA_per_cm2).
 * ALGEBRAIC[33] is a in component sarco_Cl_channel (dimensionless).
 * ALGEBRAIC[34] is J_Cl in component sarco_Cl_channel (milliV_milliM).
 * ALGEBRAIC[35] is g_Cl in component sarco_Cl_channel (milliS_per_cm2).
 * ALGEBRAIC[37] is K_R in component sarco_IR_channel (milliM).
 * ALGEBRAIC[38] is g_IR_bar in component sarco_IR_channel (milliS_per_cm2).
 * ALGEBRAIC[39] is y in component sarco_IR_channel (dimensionless).
 * ALGEBRAIC[40] is g_IR in component sarco_IR_channel (milliS_per_cm2).
 * ALGEBRAIC[1] is alpha_n in component sarco_DR_channel (per_millisecond).
 * ALGEBRAIC[14] is beta_n in component sarco_DR_channel (per_millisecond).
 * ALGEBRAIC[2] is h_K_inf in component sarco_DR_channel (dimensionless).
 * ALGEBRAIC[15] is tau_h_K in component sarco_DR_channel (millisecond).
 * STATES[8] is n in component sarco_DR_channel (dimensionless).
 * STATES[9] is h_K in component sarco_DR_channel (dimensionless).
 * ALGEBRAIC[42] is g_DR in component sarco_DR_channel (milliS_per_cm2).
 * ALGEBRAIC[3] is alpha_h in component sarco_Na_channel (per_millisecond).
 * ALGEBRAIC[16] is beta_h in component sarco_Na_channel (per_millisecond).
 * ALGEBRAIC[4] is alpha_m in component sarco_Na_channel (per_millisecond).
 * ALGEBRAIC[17] is beta_m in component sarco_Na_channel (per_millisecond).
 * ALGEBRAIC[5] is S_inf in component sarco_Na_channel (dimensionless).
 * ALGEBRAIC[18] is tau_S in component sarco_Na_channel (millisecond).
 * STATES[10] is m in component sarco_Na_channel (dimensionless).
 * STATES[11] is h in component sarco_Na_channel (dimensionless).
 * STATES[12] is S in component sarco_Na_channel (dimensionless).
 * ALGEBRAIC[45] is g_Na in component sarco_Na_channel (milliS_per_cm2).
 * ALGEBRAIC[44] is J_Na in component sarco_Na_channel (milliV_milliM).
 * ALGEBRAIC[47] is sig in component sarco_NaK_channel (dimensionless).
 * ALGEBRAIC[48] is f1 in component sarco_NaK_channel (dimensionless).
 * ALGEBRAIC[49] is I_NaK_bar in component sarco_NaK_channel (microA_per_cm2).
 * ALGEBRAIC[52] is a_t in component t_Cl_channel (dimensionless).
 * ALGEBRAIC[53] is J_Cl_t in component t_Cl_channel (milliV_milliM).
 * ALGEBRAIC[54] is g_Cl_t in component t_Cl_channel (milliS_per_cm2).
 * ALGEBRAIC[56] is K_R_t in component t_IR_channel (milliM).
 * ALGEBRAIC[57] is g_IR_bar_t in component t_IR_channel (milliS_per_cm2).
 * ALGEBRAIC[58] is y_t in component t_IR_channel (dimensionless).
 * ALGEBRAIC[59] is g_IR_t in component t_IR_channel (milliS_per_cm2).
 * ALGEBRAIC[6] is alpha_n_t in component t_DR_channel (per_millisecond).
 * ALGEBRAIC[19] is beta_n_t in component t_DR_channel (per_millisecond).
 * ALGEBRAIC[7] is h_K_inf_t in component t_DR_channel (dimensionless).
 * ALGEBRAIC[20] is tau_h_K_t in component t_DR_channel (millisecond).
 * STATES[13] is n_t in component t_DR_channel (dimensionless).
 * STATES[14] is h_K_t in component t_DR_channel (dimensionless).
 * ALGEBRAIC[61] is g_DR_t in component t_DR_channel (milliS_per_cm2).
 * ALGEBRAIC[8] is alpha_h_t in component t_Na_channel (per_millisecond).
 * ALGEBRAIC[21] is beta_h_t in component t_Na_channel (per_millisecond).
 * ALGEBRAIC[9] is alpha_m_t in component t_Na_channel (per_millisecond).
 * ALGEBRAIC[22] is beta_m_t in component t_Na_channel (per_millisecond).
 * ALGEBRAIC[10] is S_inf_t in component t_Na_channel (dimensionless).
 * ALGEBRAIC[23] is tau_S_t in component t_Na_channel (millisecond).
 * STATES[15] is m_t in component t_Na_channel (dimensionless).
 * STATES[16] is h_t in component t_Na_channel (dimensionless).
 * STATES[17] is S_t in component t_Na_channel (dimensionless).
 * ALGEBRAIC[64] is g_Na_t in component t_Na_channel (milliS_per_cm2).
 * ALGEBRAIC[63] is J_Na_t in component t_Na_channel (milliV_milliM).
 * ALGEBRAIC[66] is sig_t in component t_NaK_channel (dimensionless).
 * ALGEBRAIC[67] is f1_t in component t_NaK_channel (dimensionless).
 * ALGEBRAIC[68] is I_NaK_bar_t in component t_NaK_channel (microA_per_cm2).
 * STATES[18] is O_0 in component sternrios (dimensionless).
 * STATES[19] is O_1 in component sternrios (dimensionless).
 * STATES[20] is O_2 in component sternrios (dimensionless).
 * STATES[21] is O_3 in component sternrios (dimensionless).
 * STATES[22] is O_4 in component sternrios (dimensionless).
 * CONSTANTS[54] is k_L in component sternrios (per_millisecond).
 * CONSTANTS[55] is k_Lm in component sternrios (per_millisecond).
 * CONSTANTS[56] is f in component sternrios (dimensionless).
 * CONSTANTS[57] is alpha1 in component sternrios (per_millisecond).
 * CONSTANTS[58] is K in component sternrios (millivolt).
 * CONSTANTS[59] is Vbar in component sternrios (millivolt).
 * STATES[23] is C_0 in component sternrios (dimensionless).
 * STATES[24] is C_1 in component sternrios (dimensionless).
 * STATES[25] is C_2 in component sternrios (dimensionless).
 * STATES[26] is C_3 in component sternrios (dimensionless).
 * STATES[27] is C_4 in component sternrios (dimensionless).
 * ALGEBRAIC[11] is k_C in component sternrios (per_millisecond).
 * ALGEBRAIC[24] is k_Cm in component sternrios (per_millisecond).
 * CONSTANTS[60] is nu_SR in component razumova (micromolar_per_millisecond_micrometre3).
 * CONSTANTS[61] is K_SR in component razumova (micromolar).
 * CONSTANTS[62] is L_e in component razumova (micrometre3_per_millisecond).
 * CONSTANTS[63] is tau_R in component razumova (micrometre3_per_millisecond).
 * CONSTANTS[64] is tau_SR_R in component razumova (micrometre3_per_millisecond).
 * CONSTANTS[65] is L_x in component razumova (micrometre).
 * CONSTANTS[66] is R_R in component razumova (micrometre).
 * CONSTANTS[99] is V_o in component razumova (micrometre3).
 * CONSTANTS[101] is V_1 in component razumova (micrometre3).
 * CONSTANTS[102] is V_2 in component razumova (micrometre3).
 * CONSTANTS[100] is V_SR in component razumova (micrometre3).
 * CONSTANTS[103] is V_SR1 in component razumova (micrometre3).
 * CONSTANTS[104] is V_SR2 in component razumova (micrometre3).
 * CONSTANTS[67] is k_T_on in component razumova (per_micromolar_per_millisecond).
 * CONSTANTS[68] is k_T_off in component razumova (per_millisecond).
 * CONSTANTS[69] is T_tot in component razumova (micromolar).
 * CONSTANTS[70] is k_P_on in component razumova (per_micromolar_per_millisecond).
 * CONSTANTS[71] is k_P_off in component razumova (per_millisecond).
 * CONSTANTS[72] is P_tot in component razumova (micromolar).
 * CONSTANTS[73] is k_Mg_on in component razumova (per_micromolar_per_millisecond).
 * CONSTANTS[74] is k_Mg_off in component razumova (per_millisecond).
 * CONSTANTS[75] is k_Cs_on in component razumova (per_micromolar_per_millisecond).
 * CONSTANTS[76] is k_Cs_off in component razumova (per_millisecond).
 * CONSTANTS[77] is Cs_tot in component razumova (micromolar).
 * CONSTANTS[78] is k_CATP_on in component razumova (per_micromolar_per_millisecond).
 * CONSTANTS[79] is k_CATP_off in component razumova (per_millisecond).
 * CONSTANTS[80] is k_MATP_on in component razumova (per_micromolar_per_millisecond).
 * CONSTANTS[81] is k_MATP_off in component razumova (per_millisecond).
 * CONSTANTS[82] is tau_ATP in component razumova (micrometre3_per_millisecond).
 * CONSTANTS[83] is tau_Mg in component razumova (micrometre3_per_millisecond).
 * CONSTANTS[84] is k_0_on in component razumova (per_millisecond).
 * CONSTANTS[85] is k_0_off in component razumova (per_millisecond).
 * CONSTANTS[86] is k_Ca_on in component razumova (per_millisecond).
 * CONSTANTS[87] is k_Ca_off in component razumova (per_millisecond).
 * CONSTANTS[88] is f_o in component razumova (per_millisecond).
 * CONSTANTS[89] is f_p in component razumova (per_millisecond).
 * CONSTANTS[90] is h_o in component razumova (per_millisecond).
 * CONSTANTS[91] is h_p in component razumova (per_millisecond).
 * CONSTANTS[92] is g_o in component razumova (per_millisecond).
 * CONSTANTS[93] is b_p in component razumova (per_millisecond).
 * CONSTANTS[94] is k_p in component razumova (micrometre3_per_millisecond).
 * CONSTANTS[95] is A_p in component razumova (per_milliM3_per_millisecond).
 * CONSTANTS[96] is B_p in component razumova (per_milliM2_per_millisecond).
 * CONSTANTS[97] is PP in component razumova (milliM2).
 * ALGEBRAIC[12] is T_0 in component razumova (micromolar).
 * STATES[28] is Ca_1 in component razumova (micromolar).
 * STATES[29] is Ca_SR1 in component razumova (micromolar).
 * STATES[30] is Ca_2 in component razumova (micromolar).
 * STATES[31] is Ca_SR2 in component razumova (micromolar).
 * STATES[32] is Ca_T_2 in component razumova (micromolar).
 * STATES[33] is Ca_P1 in component razumova (micromolar).
 * STATES[34] is Ca_P2 in component razumova (micromolar).
 * STATES[35] is Mg_P1 in component razumova (micromolar).
 * STATES[36] is Mg_P2 in component razumova (micromolar).
 * STATES[37] is Ca_Cs1 in component razumova (micromolar).
 * STATES[38] is Ca_Cs2 in component razumova (micromolar).
 * STATES[39] is Ca_ATP1 in component razumova (micromolar).
 * STATES[40] is Ca_ATP2 in component razumova (micromolar).
 * STATES[41] is Mg_ATP1 in component razumova (micromolar).
 * STATES[42] is Mg_ATP2 in component razumova (micromolar).
 * STATES[43] is ATP1 in component razumova (micromolar).
 * STATES[44] is ATP2 in component razumova (micromolar).
 * STATES[45] is Mg1 in component razumova (micromolar).
 * STATES[46] is Mg2 in component razumova (micromolar).
 * STATES[47] is Ca_CaT2 in component razumova (micromolar).
 * STATES[48] is D_0 in component razumova (micromolar).
 * STATES[49] is D_1 in component razumova (micromolar).
 * STATES[50] is D_2 in component razumova (micromolar).
 * STATES[51] is A_1 in component razumova (micromolar).
 * STATES[52] is A_2 in component razumova (micromolar).
 * STATES[53] is P in component razumova (milliM).
 * STATES[54] is P_SR in component razumova (milliM).
 * STATES[55] is P_C_SR in component razumova (milliM).
 * CONSTANTS[98] is i2 in component razumova (micrometre3_per_millisecond).
 * RATES[0] is d/dt vS in component wal_environment (millivolt).
 * RATES[1] is d/dt vT in component wal_environment (millivolt).
 * RATES[3] is d/dt K_i in component wal_environment (milliM).
 * RATES[2] is d/dt K_t in component wal_environment (milliM).
 * RATES[4] is d/dt K_e in component wal_environment (milliM).
 * RATES[5] is d/dt Na_i in component wal_environment (milliM).
 * RATES[6] is d/dt Na_t in component wal_environment (milliM).
 * RATES[7] is d/dt Na_e in component wal_environment (milliM).
 * RATES[8] is d/dt n in component sarco_DR_channel (dimensionless).
 * RATES[9] is d/dt h_K in component sarco_DR_channel (dimensionless).
 * RATES[10] is d/dt m in component sarco_Na_channel (dimensionless).
 * RATES[11] is d/dt h in component sarco_Na_channel (dimensionless).
 * RATES[12] is d/dt S in component sarco_Na_channel (dimensionless).
 * RATES[13] is d/dt n_t in component t_DR_channel (dimensionless).
 * RATES[14] is d/dt h_K_t in component t_DR_channel (dimensionless).
 * RATES[15] is d/dt m_t in component t_Na_channel (dimensionless).
 * RATES[16] is d/dt h_t in component t_Na_channel (dimensionless).
 * RATES[17] is d/dt S_t in component t_Na_channel (dimensionless).
 * RATES[23] is d/dt C_0 in component sternrios (dimensionless).
 * RATES[18] is d/dt O_0 in component sternrios (dimensionless).
 * RATES[24] is d/dt C_1 in component sternrios (dimensionless).
 * RATES[19] is d/dt O_1 in component sternrios (dimensionless).
 * RATES[25] is d/dt C_2 in component sternrios (dimensionless).
 * RATES[20] is d/dt O_2 in component sternrios (dimensionless).
 * RATES[26] is d/dt C_3 in component sternrios (dimensionless).
 * RATES[21] is d/dt O_3 in component sternrios (dimensionless).
 * RATES[27] is d/dt C_4 in component sternrios (dimensionless).
 * RATES[22] is d/dt O_4 in component sternrios (dimensionless).
 * RATES[28] is d/dt Ca_1 in component razumova (micromolar).
 * RATES[29] is d/dt Ca_SR1 in component razumova (micromolar).
 * RATES[30] is d/dt Ca_2 in component razumova (micromolar).
 * RATES[31] is d/dt Ca_SR2 in component razumova (micromolar).
 * RATES[32] is d/dt Ca_T_2 in component razumova (micromolar).
 * RATES[33] is d/dt Ca_P1 in component razumova (micromolar).
 * RATES[34] is d/dt Ca_P2 in component razumova (micromolar).
 * RATES[35] is d/dt Mg_P1 in component razumova (micromolar).
 * RATES[36] is d/dt Mg_P2 in component razumova (micromolar).
 * RATES[37] is d/dt Ca_Cs1 in component razumova (micromolar).
 * RATES[38] is d/dt Ca_Cs2 in component razumova (micromolar).
 * RATES[39] is d/dt Ca_ATP1 in component razumova (micromolar).
 * RATES[40] is d/dt Ca_ATP2 in component razumova (micromolar).
 * RATES[41] is d/dt Mg_ATP1 in component razumova (micromolar).
 * RATES[42] is d/dt Mg_ATP2 in component razumova (micromolar).
 * RATES[43] is d/dt ATP1 in component razumova (micromolar).
 * RATES[44] is d/dt ATP2 in component razumova (micromolar).
 * RATES[45] is d/dt Mg1 in component razumova (micromolar).
 * RATES[46] is d/dt Mg2 in component razumova (micromolar).
 * RATES[47] is d/dt Ca_CaT2 in component razumova (micromolar).
 * RATES[48] is d/dt D_0 in component razumova (micromolar).
 * RATES[49] is d/dt D_1 in component razumova (micromolar).
 * RATES[50] is d/dt D_2 in component razumova (micromolar).
 * RATES[51] is d/dt A_1 in component razumova (micromolar).
 * RATES[52] is d/dt A_2 in component razumova (micromolar).
 * RATES[53] is d/dt P in component razumova (milliM).
 * RATES[54] is d/dt P_SR in component razumova (milliM).
 * RATES[55] is d/dt P_C_SR in component razumova (milliM).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
CONSTANTS[0] = 1.0;
CONSTANTS[1] = 4.8;
CONSTANTS[2] = 150;
CONSTANTS[3] = 0.000001;
CONSTANTS[4] = 0.0025;
CONSTANTS[5] = 0.0005;
CONSTANTS[6] = 96485;
CONSTANTS[7] = 350;
CONSTANTS[8] = 350;
CONSTANTS[9] = 0.0032;
CONSTANTS[10] = 21875;
CONSTANTS[11] = 21875;
CONSTANTS[12] = 1.02;
CONSTANTS[13] = -1.29;
CONSTANTS[14] = 0.0081;
CONSTANTS[15] = 0.288;
CONSTANTS[16] = 0.0131;
CONSTANTS[17] = 4.38;
CONSTANTS[18] = 1.38;
CONSTANTS[19] = 0.067;
CONSTANTS[20] = -46;
CONSTANTS[21] = -40;
CONSTANTS[22] = -45;
CONSTANTS[23] = 70;
CONSTANTS[24] = -78;
CONSTANTS[25] = -40;
CONSTANTS[26] = 150;
CONSTANTS[27] = 5.8;
CONSTANTS[28] = 7.5;
CONSTANTS[29] = 14.7;
CONSTANTS[30] = 9;
CONSTANTS[31] = 10;
CONSTANTS[32] = 7;
CONSTANTS[33] = 18;
CONSTANTS[34] = 40;
CONSTANTS[35] = 8314.41;
CONSTANTS[36] = 293;
CONSTANTS[37] = 19.65;
CONSTANTS[38] = 64.8;
CONSTANTS[39] = 804;
CONSTANTS[40] = 11.1;
CONSTANTS[41] = 0.4;
CONSTANTS[42] = 950;
CONSTANTS[43] = 1;
CONSTANTS[44] = 1;
CONSTANTS[45] = 13;
CONSTANTS[46] = 10;
CONSTANTS[47] = 0.000621;
CONSTANTS[48] = 90;
STATES[0] = -79.974;
STATES[1] = -80.2;
STATES[2] = 5.9;
STATES[3] = 150.9;
STATES[4] = 5.9;
STATES[5] = 12.7;
STATES[6] = 133;
STATES[7] = 133;
CONSTANTS[49] = 0.1;
CONSTANTS[50] = 1.0;
CONSTANTS[51] = 0.45;
CONSTANTS[52] = 0.1;
CONSTANTS[53] = 0.1;
STATES[8] = 0.009466;
STATES[9] = 0.9952;
STATES[10] = 0.0358;
STATES[11] = 0.4981;
STATES[12] = 0.581;
STATES[13] = 0.009466;
STATES[14] = 0.9952;
STATES[15] = 0.0358;
STATES[16] = 0.4981;
STATES[17] = 0.581;
STATES[18] = 0;
STATES[19] = 0;
STATES[20] = 0;
STATES[21] = 0;
STATES[22] = 0;
CONSTANTS[54] = 0.002;
CONSTANTS[55] = 1000;
CONSTANTS[56] = 0.2;
CONSTANTS[57] = 0.2;
CONSTANTS[58] = 4.5;
CONSTANTS[59] = -20;
STATES[23] = 1;
STATES[24] = 0;
STATES[25] = 0;
STATES[26] = 0;
STATES[27] = 0;
CONSTANTS[60] = 4.875;
CONSTANTS[61] = 1;
CONSTANTS[62] = 0.00002;
CONSTANTS[63] = 0.75;
CONSTANTS[64] = 0.75;
CONSTANTS[65] = 1.1;
CONSTANTS[66] = 0.5;
CONSTANTS[67] = 0.04425;
CONSTANTS[68] = 0.115;
CONSTANTS[69] = 140;
CONSTANTS[70] = 0.0417;
CONSTANTS[71] = 0.0005;
CONSTANTS[72] = 1500;
CONSTANTS[73] = 0.000033;
CONSTANTS[74] = 0.003;
CONSTANTS[75] = 0.000004;
CONSTANTS[76] = 0.005;
CONSTANTS[77] = 31000;
CONSTANTS[78] = 0.15;
CONSTANTS[79] = 30;
CONSTANTS[80] = 0.0015;
CONSTANTS[81] = 0.15;
CONSTANTS[82] = 0.375;
CONSTANTS[83] = 1.5;
CONSTANTS[84] = 0;
CONSTANTS[85] = 0.15;
CONSTANTS[86] = 0.15;
CONSTANTS[87] = 0.05;
CONSTANTS[88] = 1.5;
CONSTANTS[89] = 15;
CONSTANTS[90] = 0.24;
CONSTANTS[91] = 0.18;
CONSTANTS[92] = 0.12;
CONSTANTS[93] = 0.00002867;
CONSTANTS[94] = 0.00000362;
CONSTANTS[95] = 1;
CONSTANTS[96] = 0.0001;
CONSTANTS[97] = 6;
STATES[28] = 0.1;
STATES[29] = 1500;
STATES[30] = 0.1;
STATES[31] = 1500;
STATES[32] = 25;
STATES[33] = 615;
STATES[34] = 615;
STATES[35] = 811;
STATES[36] = 811;
STATES[37] = 16900;
STATES[38] = 16900;
STATES[39] = 0.4;
STATES[40] = 0.4;
STATES[41] = 7200;
STATES[42] = 7200;
STATES[43] = 799.6;
STATES[44] = 799.6;
STATES[45] = 1000;
STATES[46] = 1000;
STATES[47] = 3;
STATES[48] = 0.8;
STATES[49] = 1.2;
STATES[50] = 3;
STATES[51] = 0.3;
STATES[52] = 0.23;
STATES[53] = 0.23;
STATES[54] = 0.23;
STATES[55] = 0.23;
CONSTANTS[98] = 300;
CONSTANTS[99] =  0.950000*CONSTANTS[65]* 3.14159265358979*pow(CONSTANTS[66], 2.00000);
CONSTANTS[100] =  0.0500000*CONSTANTS[65]* 3.14159265358979*pow(CONSTANTS[66], 2.00000);
CONSTANTS[101] =  0.0100000*CONSTANTS[99];
CONSTANTS[102] =  0.990000*CONSTANTS[99];
CONSTANTS[103] =  0.0100000*CONSTANTS[100];
CONSTANTS[104] =  0.990000*CONSTANTS[100];
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[28] = (((( ( CONSTANTS[98]*(STATES[18]+STATES[19]+STATES[20]+STATES[21]+STATES[22]))*((STATES[29] - STATES[28])/CONSTANTS[101]) -  CONSTANTS[60]*((STATES[28]/(STATES[28]+CONSTANTS[61]))/CONSTANTS[101]))+ CONSTANTS[62]*((STATES[29] - STATES[28])/CONSTANTS[101]))+ - CONSTANTS[63]*((STATES[28] - STATES[30])/CONSTANTS[101]))+- ( ( CONSTANTS[70]*STATES[28])*((CONSTANTS[72]+- STATES[33])+- STATES[35])+ - CONSTANTS[71]*STATES[33]))+- ( ( CONSTANTS[78]*STATES[28])*STATES[43]+ - CONSTANTS[79]*STATES[39]);
RATES[29] = ((( - ( CONSTANTS[98]*(STATES[18]+STATES[19]+STATES[20]+STATES[21]+STATES[22]))*((STATES[29] - STATES[28])/CONSTANTS[103])+ CONSTANTS[60]*((STATES[28]/(STATES[28]+CONSTANTS[61]))/CONSTANTS[103]))+ - CONSTANTS[62]*((STATES[29] - STATES[28])/CONSTANTS[103]))+ - CONSTANTS[64]*((STATES[29] - STATES[31])/CONSTANTS[103]))+- ( ( CONSTANTS[75]*STATES[29])*(CONSTANTS[77] - STATES[37])+ - CONSTANTS[76]*STATES[37]);
RATES[31] = ((( CONSTANTS[60]*((STATES[30]/(STATES[30]+CONSTANTS[61]))/CONSTANTS[104])+ - CONSTANTS[62]*((STATES[31]+- STATES[30])/CONSTANTS[104]))+ CONSTANTS[64]*((STATES[29]+- STATES[31])/CONSTANTS[104]))+- ( ( CONSTANTS[75]*STATES[31])*(CONSTANTS[77]+- STATES[38])+ - CONSTANTS[76]*STATES[38])) -  (1000.00/1.00000)*( CONSTANTS[95]*( STATES[54]*(0.00100000/1.00000)*STATES[31] - CONSTANTS[97])*( STATES[54]*(0.00100000/1.00000)*STATES[31] - CONSTANTS[97]>0.00000 ? 1.00000 : 0.00000)*(0.00100000/1.00000)*STATES[54]*STATES[31] -  CONSTANTS[96]*STATES[55]*(CONSTANTS[97] -  STATES[54]*(0.00100000/1.00000)*STATES[31])*(CONSTANTS[97] -  STATES[54]*(0.00100000/1.00000)*STATES[31]>0.00000 ? 1.00000 : 0.00000));
RATES[33] =  ( CONSTANTS[70]*STATES[28])*((CONSTANTS[72]+- STATES[33])+- STATES[35])+ - CONSTANTS[71]*STATES[33];
RATES[34] =  ( CONSTANTS[70]*STATES[30])*((CONSTANTS[72]+- STATES[34])+- STATES[36])+ - CONSTANTS[71]*STATES[34];
RATES[35] =  ( CONSTANTS[73]*(CONSTANTS[72]+- STATES[33]+- STATES[35]))*STATES[45]+ - CONSTANTS[74]*STATES[35];
RATES[36] =  ( CONSTANTS[73]*(CONSTANTS[72]+- STATES[34]+- STATES[36]))*STATES[46]+ - CONSTANTS[74]*STATES[36];
RATES[37] =  ( CONSTANTS[75]*STATES[29])*(CONSTANTS[77]+- STATES[37])+ - CONSTANTS[76]*STATES[37];
RATES[38] =  ( CONSTANTS[75]*STATES[31])*(CONSTANTS[77]+- STATES[38])+ - CONSTANTS[76]*STATES[38];
RATES[39] = ( ( CONSTANTS[78]*STATES[28])*STATES[43]+ - CONSTANTS[79]*STATES[39])+ - CONSTANTS[82]*((STATES[39]+- STATES[40])/CONSTANTS[101]);
RATES[40] = ( ( CONSTANTS[78]*STATES[30])*STATES[44]+ - CONSTANTS[79]*STATES[40])+ CONSTANTS[82]*((STATES[39]+- STATES[40])/CONSTANTS[102]);
RATES[41] = ( ( CONSTANTS[80]*STATES[45])*STATES[43]+ - CONSTANTS[81]*STATES[41])+ - CONSTANTS[82]*((STATES[41]+- STATES[42])/CONSTANTS[101]);
RATES[42] = ( ( CONSTANTS[80]*STATES[46])*STATES[44]+ - CONSTANTS[81]*STATES[42])+ CONSTANTS[82]*((STATES[41]+- STATES[42])/CONSTANTS[102]);
RATES[43] = (- ( ( CONSTANTS[78]*STATES[28])*STATES[43]+ - CONSTANTS[79]*STATES[39])+- ( ( CONSTANTS[80]*STATES[45])*STATES[43]+ - CONSTANTS[81]*STATES[41]))+ - CONSTANTS[82]*((STATES[43]+- STATES[44])/CONSTANTS[101]);
RATES[44] = (- ( ( CONSTANTS[78]*STATES[30])*STATES[44]+ - CONSTANTS[79]*STATES[40])+- ( ( CONSTANTS[80]*STATES[46])*STATES[44]+ - CONSTANTS[81]*STATES[42]))+ CONSTANTS[82]*((STATES[43]+- STATES[44])/CONSTANTS[102]);
RATES[45] = (- ( ( CONSTANTS[73]*(CONSTANTS[72]+- STATES[33]+- STATES[35]))*STATES[45]+ - CONSTANTS[74]*STATES[35])+- ( ( CONSTANTS[80]*STATES[45])*STATES[43]+ - CONSTANTS[81]*STATES[41]))+ - CONSTANTS[83]*((STATES[45]+- STATES[46])/CONSTANTS[101]);
RATES[46] = (- ( ( CONSTANTS[73]*(CONSTANTS[72]+- STATES[34]+- STATES[36]))*STATES[46]+ - CONSTANTS[74]*STATES[36])+- ( ( CONSTANTS[80]*STATES[46])*STATES[44]+ - CONSTANTS[81]*STATES[42]))+ CONSTANTS[83]*((STATES[45]+- STATES[46])/CONSTANTS[102]);
RATES[47] = (( ( CONSTANTS[67]*STATES[30])*STATES[32]+ - CONSTANTS[68]*STATES[47])+ - CONSTANTS[86]*STATES[47])+ CONSTANTS[87]*STATES[50];
RATES[49] = (((( CONSTANTS[67]*STATES[30]*STATES[48]+ - CONSTANTS[68]*STATES[49])+ CONSTANTS[84]*STATES[32])+ - CONSTANTS[85]*STATES[49])+ ( - CONSTANTS[67]*STATES[30])*STATES[49])+ CONSTANTS[68]*STATES[50];
RATES[50] = ((((( CONSTANTS[67]*STATES[30]*STATES[49]+ - CONSTANTS[68]*STATES[50])+ CONSTANTS[86]*STATES[47])+ - CONSTANTS[87]*STATES[50])+ - CONSTANTS[88]*STATES[50])+ CONSTANTS[89]*STATES[51])+ CONSTANTS[92]*STATES[52];
RATES[51] = (( CONSTANTS[88]*STATES[50]+ - CONSTANTS[89]*STATES[51])+ CONSTANTS[91]*STATES[52])+ - CONSTANTS[90]*STATES[51];
RATES[52] = ( - CONSTANTS[91]*STATES[52]+ CONSTANTS[90]*STATES[51])+ - CONSTANTS[92]*STATES[52];
RATES[53] =  (0.00100000/1.00000)*( CONSTANTS[90]*STATES[51] -  CONSTANTS[91]*STATES[52])+ -1.00000*CONSTANTS[93]*STATES[53]+ -1.00000*CONSTANTS[94]*((STATES[53] - STATES[54])/CONSTANTS[102]);
RATES[54] =  CONSTANTS[94]*((STATES[53] - STATES[54])/CONSTANTS[104]) -  1.00000*( CONSTANTS[95]*( STATES[54]*(0.00100000/1.00000)*STATES[31] - CONSTANTS[97])*( STATES[54]*(0.00100000/1.00000)*STATES[31] - CONSTANTS[97]>0.00000 ? 1.00000 : 0.00000)*(0.00100000/1.00000)*STATES[54]*STATES[31] -  CONSTANTS[96]*STATES[55]*(CONSTANTS[97] -  STATES[54]*(0.00100000/1.00000)*STATES[31])*(CONSTANTS[97] -  STATES[54]*(0.00100000/1.00000)*STATES[31]>0.00000 ? 1.00000 : 0.00000));
RATES[55] =  1.00000*( CONSTANTS[95]*( STATES[54]*(0.00100000/1.00000)*STATES[31] - CONSTANTS[97])*( STATES[54]*(0.00100000/1.00000)*STATES[31] - CONSTANTS[97]>0.00000 ? 1.00000 : 0.00000)*(0.00100000/1.00000)*STATES[54]*STATES[31] -  CONSTANTS[96]*STATES[55]*(CONSTANTS[97] -  STATES[54]*(0.00100000/1.00000)*STATES[31])*(CONSTANTS[97] -  STATES[54]*(0.00100000/1.00000)*STATES[31]>0.00000 ? 1.00000 : 0.00000));
ALGEBRAIC[12] = CONSTANTS[69]+- STATES[32]+- STATES[47]+- STATES[48]+- STATES[49]+- STATES[50]+- STATES[51]+- STATES[52];
RATES[30] = (((( - CONSTANTS[60]*((STATES[30]/(STATES[30]+CONSTANTS[61]))/CONSTANTS[102])+ CONSTANTS[62]*((STATES[31]+- STATES[30])/CONSTANTS[102]))+ CONSTANTS[63]*((STATES[28] - STATES[30])/CONSTANTS[102]))+- ((((((( CONSTANTS[67]*STATES[30]*ALGEBRAIC[12]+ - CONSTANTS[68]*STATES[32])+ CONSTANTS[67]*STATES[30]*STATES[32])+ - CONSTANTS[68]*STATES[47])+ CONSTANTS[67]*STATES[30]*STATES[48])+ - CONSTANTS[68]*STATES[49])+ CONSTANTS[67]*STATES[30]*STATES[49])+ - CONSTANTS[68]*STATES[50]))+- ( ( CONSTANTS[70]*STATES[30])*(CONSTANTS[72]+- STATES[34]+- STATES[36])+ - CONSTANTS[71]*STATES[34]))+- ( ( CONSTANTS[78]*STATES[30])*STATES[44]+ - CONSTANTS[79]*STATES[40]);
RATES[32] = (((( ( CONSTANTS[67]*STATES[30])*ALGEBRAIC[12]+ - CONSTANTS[68]*STATES[32])+ ( - CONSTANTS[67]*STATES[30])*STATES[32])+ CONSTANTS[68]*STATES[47])+ - CONSTANTS[84]*STATES[32])+ CONSTANTS[85]*STATES[49];
RATES[48] = (( ( - CONSTANTS[67]*STATES[30])*STATES[48]+ CONSTANTS[68]*STATES[49])+ CONSTANTS[84]*ALGEBRAIC[12])+ - CONSTANTS[85]*STATES[48];
ALGEBRAIC[1] =  CONSTANTS[16]*((STATES[0] - CONSTANTS[21])/(1.00000 - exp(- ((STATES[0] - CONSTANTS[21])/CONSTANTS[32]))));
ALGEBRAIC[14] =  CONSTANTS[19]*exp(- ((STATES[0] - CONSTANTS[21])/CONSTANTS[34]));
RATES[8] =  ALGEBRAIC[1]*(1.00000 - STATES[8]) -  ALGEBRAIC[14]*STATES[8];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0] - CONSTANTS[25])/CONSTANTS[28]));
ALGEBRAIC[15] =  1000.00*exp(- ((STATES[0]+40.0000)/25.7500));
RATES[9] = (ALGEBRAIC[2] - STATES[9])/ALGEBRAIC[15];
ALGEBRAIC[4] =  CONSTANTS[15]*((STATES[0] - CONSTANTS[20])/(1.00000 - exp(- ((STATES[0] - CONSTANTS[20])/CONSTANTS[31]))));
ALGEBRAIC[17] =  CONSTANTS[18]*exp(- ((STATES[0] - CONSTANTS[20])/CONSTANTS[33]));
RATES[10] =  ALGEBRAIC[4]*(1.00000 - STATES[10]) -  ALGEBRAIC[17]*STATES[10];
ALGEBRAIC[3] =  CONSTANTS[14]*exp(- ((STATES[0] - CONSTANTS[22])/CONSTANTS[29]));
ALGEBRAIC[16] = CONSTANTS[17]/(1.00000+exp(- ((STATES[0] - CONSTANTS[22])/CONSTANTS[30])));
RATES[11] =  ALGEBRAIC[3]*(1.00000 - STATES[11]) -  ALGEBRAIC[16]*STATES[11];
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0] - CONSTANTS[24])/CONSTANTS[27]));
ALGEBRAIC[18] = 8571.00/(0.200000+ 5.65000*pow((STATES[0]+CONSTANTS[48])/100.000, 2.00000));
RATES[12] = (ALGEBRAIC[5] - STATES[12])/ALGEBRAIC[18];
ALGEBRAIC[6] =  CONSTANTS[16]*((STATES[1] - CONSTANTS[21])/(1.00000 - exp(- ((STATES[1] - CONSTANTS[21])/CONSTANTS[32]))));
ALGEBRAIC[19] =  CONSTANTS[19]*exp(- ((STATES[1] - CONSTANTS[21])/CONSTANTS[34]));
RATES[13] =  ALGEBRAIC[6]*(1.00000 - STATES[13]) -  ALGEBRAIC[19]*STATES[13];
ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[1] - CONSTANTS[25])/CONSTANTS[28]));
ALGEBRAIC[20] =  1.00000*exp(- ((STATES[1]+40.0000)/25.7500));
RATES[14] = (ALGEBRAIC[7] - STATES[14])/ALGEBRAIC[20];
ALGEBRAIC[9] =  CONSTANTS[15]*((STATES[1] - CONSTANTS[20])/(1.00000 - exp(- ((STATES[1] - CONSTANTS[20])/CONSTANTS[31]))));
ALGEBRAIC[22] =  CONSTANTS[18]*exp(- ((STATES[1] - CONSTANTS[20])/CONSTANTS[33]));
RATES[15] =  ALGEBRAIC[9]*(1.00000 - STATES[15]) -  ALGEBRAIC[22]*STATES[15];
ALGEBRAIC[8] =  CONSTANTS[14]*exp(- ((STATES[1] - CONSTANTS[22])/CONSTANTS[29]));
ALGEBRAIC[21] = CONSTANTS[17]/(1.00000+exp(- ((STATES[1] - CONSTANTS[22])/CONSTANTS[30])));
RATES[16] =  ALGEBRAIC[8]*(1.00000 - STATES[16]) -  ALGEBRAIC[21]*STATES[16];
ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[1] - CONSTANTS[24])/CONSTANTS[27]));
ALGEBRAIC[23] = 8571.00/(0.200000+ 5.65000*pow((STATES[1]+CONSTANTS[48])/100.000, 2.00000));
RATES[17] = (ALGEBRAIC[10] - STATES[17])/ALGEBRAIC[23];
ALGEBRAIC[11] =  0.500000*CONSTANTS[57]*exp((STATES[1] - CONSTANTS[59])/( 8.00000*CONSTANTS[58]));
ALGEBRAIC[24] =  0.500000*CONSTANTS[57]*exp((CONSTANTS[59] - STATES[1])/( 8.00000*CONSTANTS[58]));
RATES[23] =  - CONSTANTS[54]*STATES[23]+ CONSTANTS[55]*STATES[18]+ -4.00000*ALGEBRAIC[11]*STATES[23]+ ALGEBRAIC[24]*STATES[24];
RATES[18] =  CONSTANTS[54]*STATES[23]+ - CONSTANTS[55]*STATES[18]+( -4.00000*ALGEBRAIC[11]*STATES[18])/CONSTANTS[56]+ CONSTANTS[56]*ALGEBRAIC[24]*STATES[19];
RATES[24] =  4.00000*ALGEBRAIC[11]*STATES[23]+ - ALGEBRAIC[24]*STATES[24]+( - CONSTANTS[54]*STATES[24])/CONSTANTS[56]+ CONSTANTS[56]*CONSTANTS[55]*STATES[19]+ -3.00000*ALGEBRAIC[11]*STATES[24]+ 2.00000*ALGEBRAIC[24]*STATES[25];
RATES[19] = ( CONSTANTS[54]*STATES[24])/CONSTANTS[56]+ - CONSTANTS[55]*CONSTANTS[56]*STATES[19]+( 4.00000*ALGEBRAIC[11]*STATES[18])/CONSTANTS[56]+ - CONSTANTS[56]*ALGEBRAIC[24]*STATES[19]+( -3.00000*ALGEBRAIC[11]*STATES[19])/CONSTANTS[56]+ 2.00000*CONSTANTS[56]*ALGEBRAIC[24]*STATES[20];
RATES[25] =  3.00000*ALGEBRAIC[11]*STATES[24]+ -2.00000*ALGEBRAIC[24]*STATES[25]+( - CONSTANTS[54]*STATES[25])/pow(CONSTANTS[56], 2.00000)+ pow(CONSTANTS[56], 2.00000)*CONSTANTS[55]*STATES[20]+ -2.00000*ALGEBRAIC[11]*STATES[25]+ 3.00000*ALGEBRAIC[24]*STATES[26];
RATES[20] = ( 3.00000*ALGEBRAIC[11]*STATES[19])/CONSTANTS[56]+ -2.00000*CONSTANTS[56]*ALGEBRAIC[24]*STATES[20]+( CONSTANTS[54]*STATES[25])/pow(CONSTANTS[56], 2.00000)+ - CONSTANTS[55]*pow(CONSTANTS[56], 2.00000)*STATES[20]+( -2.00000*ALGEBRAIC[11]*STATES[20])/CONSTANTS[56]+ 3.00000*CONSTANTS[56]*ALGEBRAIC[24]*STATES[21];
RATES[26] =  2.00000*ALGEBRAIC[11]*STATES[25]+ -3.00000*ALGEBRAIC[24]*STATES[26]+( - CONSTANTS[54]*STATES[26])/pow(CONSTANTS[56], 3.00000)+ CONSTANTS[55]*pow(CONSTANTS[56], 3.00000)*STATES[21]+ - ALGEBRAIC[11]*STATES[26]+ 4.00000*ALGEBRAIC[24]*STATES[27];
RATES[21] = ( CONSTANTS[54]*STATES[26])/pow(CONSTANTS[56], 3.00000)+ - CONSTANTS[55]*pow(CONSTANTS[56], 3.00000)*STATES[21]+( 2.00000*ALGEBRAIC[11]*STATES[20])/CONSTANTS[56]+ -3.00000*ALGEBRAIC[24]*CONSTANTS[56]*STATES[21]+( - ALGEBRAIC[11]*STATES[21])/CONSTANTS[56]+ 4.00000*CONSTANTS[56]*ALGEBRAIC[24]*STATES[22];
RATES[27] =  ALGEBRAIC[11]*STATES[26]+ -4.00000*ALGEBRAIC[24]*STATES[27]+( - CONSTANTS[54]*STATES[27])/pow(CONSTANTS[56], 4.00000)+ CONSTANTS[55]*pow(CONSTANTS[56], 4.00000)*STATES[22];
RATES[22] = ( ALGEBRAIC[11]*STATES[21])/CONSTANTS[56]+ -4.00000*CONSTANTS[56]*ALGEBRAIC[24]*STATES[22]+( CONSTANTS[54]*STATES[27])/pow(CONSTANTS[56], 4.00000)+ - CONSTANTS[55]*pow(CONSTANTS[56], 4.00000)*STATES[22];
ALGEBRAIC[30] =  STATES[0]*((STATES[3] -  STATES[4]*exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[13] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(STATES[4]/STATES[3]);
ALGEBRAIC[37] =  STATES[4]*exp( ( - CONSTANTS[41]*ALGEBRAIC[13])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[38] =  CONSTANTS[40]*(pow(ALGEBRAIC[37], 2.00000)/(CONSTANTS[42]+pow(ALGEBRAIC[37], 2.00000)));
ALGEBRAIC[39] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(ALGEBRAIC[37], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*STATES[0]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
ALGEBRAIC[40] =  ALGEBRAIC[38]*ALGEBRAIC[39];
ALGEBRAIC[41] =  ALGEBRAIC[40]*(ALGEBRAIC[30]>0.00000 ? 1.00000 : 0.00000)*(ALGEBRAIC[30]/50.0000);
ALGEBRAIC[42] =  ( CONSTANTS[38]*pow(STATES[8], 4.00000))*STATES[9];
ALGEBRAIC[43] =  ALGEBRAIC[42]*(ALGEBRAIC[30]/50.0000);
ALGEBRAIC[47] =  (1.00000/7.00000)*(exp(STATES[7]/67.3000) - 1.00000);
ALGEBRAIC[48] = pow(1.00000+ 0.120000*exp( -0.100000*STATES[0]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*ALGEBRAIC[47]*exp(- ( STATES[0]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
ALGEBRAIC[49] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/STATES[4], 2.00000)*pow(1.00000+CONSTANTS[45]/STATES[5], 3.00000)));
ALGEBRAIC[50] =  ALGEBRAIC[49]*ALGEBRAIC[48];
RATES[4] = (ALGEBRAIC[41]+ALGEBRAIC[43]+CONSTANTS[12]+ - 2.00000*ALGEBRAIC[50])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[5])+(STATES[2] - STATES[4])/CONSTANTS[10];
ALGEBRAIC[45] =  ( ( CONSTANTS[39]*pow(STATES[10], 3.00000))*STATES[11])*STATES[12];
ALGEBRAIC[44] =  STATES[0]*((STATES[5] -  STATES[7]*exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[46] =  ALGEBRAIC[45]*(ALGEBRAIC[44]/75.0000);
RATES[7] = (ALGEBRAIC[46]+CONSTANTS[13]+ 3.00000*ALGEBRAIC[50])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[5])+(STATES[6] - STATES[7])/CONSTANTS[11];
ALGEBRAIC[0] =  (1000.00/1.00000)*((STATES[0] - STATES[1])/CONSTANTS[2]);
ALGEBRAIC[26] = 156.500/(5.00000+exp(( - CONSTANTS[6]*ALGEBRAIC[13])/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[27] = 156.500 -  5.00000*ALGEBRAIC[26];
ALGEBRAIC[34] =  STATES[0]*((ALGEBRAIC[26] -  ALGEBRAIC[27]*exp(( CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[33] = 1.00000/(1.00000+exp((STATES[0] - CONSTANTS[23])/CONSTANTS[26]));
ALGEBRAIC[35] =  CONSTANTS[37]*pow(ALGEBRAIC[33], 4.00000);
ALGEBRAIC[36] =  ALGEBRAIC[35]*(ALGEBRAIC[34]/45.0000);
ALGEBRAIC[32] = (VOI>=0.00000&&VOI<0.500000 ? 150.000 : VOI>=50.0000&&VOI<50.5000 ? 150.000 : VOI>=100.000&&VOI<100.500 ? 150.000 : VOI>=150.000&&VOI<150.500 ? 150.000 : VOI>=200.000&&VOI<200.500 ? 150.000 : VOI>=250.000&&VOI<250.500 ? 150.000 : VOI>=300.000&&VOI<300.500 ? 150.000 : VOI>=350.000&&VOI<350.500 ? 150.000 : VOI>=400.000&&VOI<400.500 ? 150.000 : 0.00000);
ALGEBRAIC[51] = ALGEBRAIC[36]+ALGEBRAIC[41]+ALGEBRAIC[43]+ALGEBRAIC[46]+ALGEBRAIC[50]+- ALGEBRAIC[32];
RATES[0] = - ((ALGEBRAIC[51]+ALGEBRAIC[0])/CONSTANTS[0]);
ALGEBRAIC[31] =  STATES[1]*((STATES[3] -  STATES[2]*exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[25] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(STATES[2]/STATES[3]);
ALGEBRAIC[56] =  STATES[2]*exp( ( - CONSTANTS[41]*ALGEBRAIC[25])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[57] =  CONSTANTS[40]*(pow(ALGEBRAIC[56], 2.00000)/(CONSTANTS[42]+pow(ALGEBRAIC[56], 2.00000)));
ALGEBRAIC[58] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(ALGEBRAIC[56], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*STATES[1]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
ALGEBRAIC[59] =  ALGEBRAIC[57]*ALGEBRAIC[58];
ALGEBRAIC[60] =  CONSTANTS[50]*ALGEBRAIC[59]*(ALGEBRAIC[31]/50.0000);
ALGEBRAIC[61] =  ( CONSTANTS[38]*pow(STATES[13], 4.00000))*STATES[14];
ALGEBRAIC[62] =  CONSTANTS[51]*ALGEBRAIC[61]*(ALGEBRAIC[31]/50.0000);
ALGEBRAIC[66] =  (1.00000/7.00000)*(exp(STATES[6]/67.3000) - 1.00000);
ALGEBRAIC[67] = pow(1.00000+ 0.120000*exp( -0.100000*STATES[1]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*ALGEBRAIC[66]*exp(- ( STATES[1]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
ALGEBRAIC[68] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/STATES[2], 2.00000)*pow(1.00000+CONSTANTS[45]/STATES[5], 3.00000)));
ALGEBRAIC[69] =  CONSTANTS[53]*ALGEBRAIC[68]*ALGEBRAIC[67];
RATES[3] =  - CONSTANTS[9]*((ALGEBRAIC[60]+ALGEBRAIC[62]+CONSTANTS[12]+ - 2.00000*ALGEBRAIC[69])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3])) - (ALGEBRAIC[41]+ALGEBRAIC[43]+CONSTANTS[12]+ -2.00000*ALGEBRAIC[50])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[4]);
RATES[2] = (ALGEBRAIC[60]+ALGEBRAIC[62]+CONSTANTS[12]+ - 2.00000*ALGEBRAIC[69])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3]) - (STATES[2] - STATES[4])/CONSTANTS[7];
ALGEBRAIC[64] =  ( ( CONSTANTS[39]*pow(STATES[15], 3.00000))*STATES[16])*STATES[17];
ALGEBRAIC[63] =  STATES[1]*((STATES[5] -  STATES[6]*exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[65] =  CONSTANTS[52]*ALGEBRAIC[64]*(ALGEBRAIC[63]/75.0000);
RATES[5] =  - CONSTANTS[9]*((ALGEBRAIC[65]+CONSTANTS[13]+ 3.00000*ALGEBRAIC[69])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3])) - (ALGEBRAIC[46]+CONSTANTS[13]+ 3.00000*ALGEBRAIC[50])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[4]);
RATES[6] = (ALGEBRAIC[65]+CONSTANTS[13]+ 3.00000*ALGEBRAIC[69])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3]) - (STATES[6] - STATES[7])/CONSTANTS[8];
ALGEBRAIC[28] = 156.500/(5.00000+exp(( - CONSTANTS[6]*ALGEBRAIC[25])/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[29] = 156.500 -  5.00000*ALGEBRAIC[28];
ALGEBRAIC[53] =  STATES[1]*((ALGEBRAIC[28] -  ALGEBRAIC[29]*exp(( CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[52] = 1.00000/(1.00000+exp((STATES[1] - CONSTANTS[23])/CONSTANTS[26]));
ALGEBRAIC[54] =  CONSTANTS[37]*pow(ALGEBRAIC[52], 4.00000);
ALGEBRAIC[55] =  CONSTANTS[49]*ALGEBRAIC[54]*(ALGEBRAIC[53]/45.0000);
ALGEBRAIC[70] = ALGEBRAIC[55]+ALGEBRAIC[60]+ALGEBRAIC[62]+ALGEBRAIC[65]+ALGEBRAIC[69];
RATES[1] = - ((ALGEBRAIC[70] - ALGEBRAIC[0]/CONSTANTS[1])/CONSTANTS[0]);
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[12] = CONSTANTS[69]+- STATES[32]+- STATES[47]+- STATES[48]+- STATES[49]+- STATES[50]+- STATES[51]+- STATES[52];
ALGEBRAIC[1] =  CONSTANTS[16]*((STATES[0] - CONSTANTS[21])/(1.00000 - exp(- ((STATES[0] - CONSTANTS[21])/CONSTANTS[32]))));
ALGEBRAIC[14] =  CONSTANTS[19]*exp(- ((STATES[0] - CONSTANTS[21])/CONSTANTS[34]));
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0] - CONSTANTS[25])/CONSTANTS[28]));
ALGEBRAIC[15] =  1000.00*exp(- ((STATES[0]+40.0000)/25.7500));
ALGEBRAIC[4] =  CONSTANTS[15]*((STATES[0] - CONSTANTS[20])/(1.00000 - exp(- ((STATES[0] - CONSTANTS[20])/CONSTANTS[31]))));
ALGEBRAIC[17] =  CONSTANTS[18]*exp(- ((STATES[0] - CONSTANTS[20])/CONSTANTS[33]));
ALGEBRAIC[3] =  CONSTANTS[14]*exp(- ((STATES[0] - CONSTANTS[22])/CONSTANTS[29]));
ALGEBRAIC[16] = CONSTANTS[17]/(1.00000+exp(- ((STATES[0] - CONSTANTS[22])/CONSTANTS[30])));
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0] - CONSTANTS[24])/CONSTANTS[27]));
ALGEBRAIC[18] = 8571.00/(0.200000+ 5.65000*pow((STATES[0]+CONSTANTS[48])/100.000, 2.00000));
ALGEBRAIC[6] =  CONSTANTS[16]*((STATES[1] - CONSTANTS[21])/(1.00000 - exp(- ((STATES[1] - CONSTANTS[21])/CONSTANTS[32]))));
ALGEBRAIC[19] =  CONSTANTS[19]*exp(- ((STATES[1] - CONSTANTS[21])/CONSTANTS[34]));
ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[1] - CONSTANTS[25])/CONSTANTS[28]));
ALGEBRAIC[20] =  1.00000*exp(- ((STATES[1]+40.0000)/25.7500));
ALGEBRAIC[9] =  CONSTANTS[15]*((STATES[1] - CONSTANTS[20])/(1.00000 - exp(- ((STATES[1] - CONSTANTS[20])/CONSTANTS[31]))));
ALGEBRAIC[22] =  CONSTANTS[18]*exp(- ((STATES[1] - CONSTANTS[20])/CONSTANTS[33]));
ALGEBRAIC[8] =  CONSTANTS[14]*exp(- ((STATES[1] - CONSTANTS[22])/CONSTANTS[29]));
ALGEBRAIC[21] = CONSTANTS[17]/(1.00000+exp(- ((STATES[1] - CONSTANTS[22])/CONSTANTS[30])));
ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[1] - CONSTANTS[24])/CONSTANTS[27]));
ALGEBRAIC[23] = 8571.00/(0.200000+ 5.65000*pow((STATES[1]+CONSTANTS[48])/100.000, 2.00000));
ALGEBRAIC[11] =  0.500000*CONSTANTS[57]*exp((STATES[1] - CONSTANTS[59])/( 8.00000*CONSTANTS[58]));
ALGEBRAIC[24] =  0.500000*CONSTANTS[57]*exp((CONSTANTS[59] - STATES[1])/( 8.00000*CONSTANTS[58]));
ALGEBRAIC[30] =  STATES[0]*((STATES[3] -  STATES[4]*exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[13] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(STATES[4]/STATES[3]);
ALGEBRAIC[37] =  STATES[4]*exp( ( - CONSTANTS[41]*ALGEBRAIC[13])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[38] =  CONSTANTS[40]*(pow(ALGEBRAIC[37], 2.00000)/(CONSTANTS[42]+pow(ALGEBRAIC[37], 2.00000)));
ALGEBRAIC[39] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(ALGEBRAIC[37], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*STATES[0]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
ALGEBRAIC[40] =  ALGEBRAIC[38]*ALGEBRAIC[39];
ALGEBRAIC[41] =  ALGEBRAIC[40]*(ALGEBRAIC[30]>0.00000 ? 1.00000 : 0.00000)*(ALGEBRAIC[30]/50.0000);
ALGEBRAIC[42] =  ( CONSTANTS[38]*pow(STATES[8], 4.00000))*STATES[9];
ALGEBRAIC[43] =  ALGEBRAIC[42]*(ALGEBRAIC[30]/50.0000);
ALGEBRAIC[47] =  (1.00000/7.00000)*(exp(STATES[7]/67.3000) - 1.00000);
ALGEBRAIC[48] = pow(1.00000+ 0.120000*exp( -0.100000*STATES[0]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*ALGEBRAIC[47]*exp(- ( STATES[0]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
ALGEBRAIC[49] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/STATES[4], 2.00000)*pow(1.00000+CONSTANTS[45]/STATES[5], 3.00000)));
ALGEBRAIC[50] =  ALGEBRAIC[49]*ALGEBRAIC[48];
ALGEBRAIC[45] =  ( ( CONSTANTS[39]*pow(STATES[10], 3.00000))*STATES[11])*STATES[12];
ALGEBRAIC[44] =  STATES[0]*((STATES[5] -  STATES[7]*exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[46] =  ALGEBRAIC[45]*(ALGEBRAIC[44]/75.0000);
ALGEBRAIC[0] =  (1000.00/1.00000)*((STATES[0] - STATES[1])/CONSTANTS[2]);
ALGEBRAIC[26] = 156.500/(5.00000+exp(( - CONSTANTS[6]*ALGEBRAIC[13])/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[27] = 156.500 -  5.00000*ALGEBRAIC[26];
ALGEBRAIC[34] =  STATES[0]*((ALGEBRAIC[26] -  ALGEBRAIC[27]*exp(( CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*STATES[0])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[33] = 1.00000/(1.00000+exp((STATES[0] - CONSTANTS[23])/CONSTANTS[26]));
ALGEBRAIC[35] =  CONSTANTS[37]*pow(ALGEBRAIC[33], 4.00000);
ALGEBRAIC[36] =  ALGEBRAIC[35]*(ALGEBRAIC[34]/45.0000);
ALGEBRAIC[32] = (VOI>=0.00000&&VOI<0.500000 ? 150.000 : VOI>=50.0000&&VOI<50.5000 ? 150.000 : VOI>=100.000&&VOI<100.500 ? 150.000 : VOI>=150.000&&VOI<150.500 ? 150.000 : VOI>=200.000&&VOI<200.500 ? 150.000 : VOI>=250.000&&VOI<250.500 ? 150.000 : VOI>=300.000&&VOI<300.500 ? 150.000 : VOI>=350.000&&VOI<350.500 ? 150.000 : VOI>=400.000&&VOI<400.500 ? 150.000 : 0.00000);
ALGEBRAIC[51] = ALGEBRAIC[36]+ALGEBRAIC[41]+ALGEBRAIC[43]+ALGEBRAIC[46]+ALGEBRAIC[50]+- ALGEBRAIC[32];
ALGEBRAIC[31] =  STATES[1]*((STATES[3] -  STATES[2]*exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[25] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(STATES[2]/STATES[3]);
ALGEBRAIC[56] =  STATES[2]*exp( ( - CONSTANTS[41]*ALGEBRAIC[25])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[57] =  CONSTANTS[40]*(pow(ALGEBRAIC[56], 2.00000)/(CONSTANTS[42]+pow(ALGEBRAIC[56], 2.00000)));
ALGEBRAIC[58] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(ALGEBRAIC[56], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*STATES[1]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
ALGEBRAIC[59] =  ALGEBRAIC[57]*ALGEBRAIC[58];
ALGEBRAIC[60] =  CONSTANTS[50]*ALGEBRAIC[59]*(ALGEBRAIC[31]/50.0000);
ALGEBRAIC[61] =  ( CONSTANTS[38]*pow(STATES[13], 4.00000))*STATES[14];
ALGEBRAIC[62] =  CONSTANTS[51]*ALGEBRAIC[61]*(ALGEBRAIC[31]/50.0000);
ALGEBRAIC[66] =  (1.00000/7.00000)*(exp(STATES[6]/67.3000) - 1.00000);
ALGEBRAIC[67] = pow(1.00000+ 0.120000*exp( -0.100000*STATES[1]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*ALGEBRAIC[66]*exp(- ( STATES[1]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
ALGEBRAIC[68] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/STATES[2], 2.00000)*pow(1.00000+CONSTANTS[45]/STATES[5], 3.00000)));
ALGEBRAIC[69] =  CONSTANTS[53]*ALGEBRAIC[68]*ALGEBRAIC[67];
ALGEBRAIC[64] =  ( ( CONSTANTS[39]*pow(STATES[15], 3.00000))*STATES[16])*STATES[17];
ALGEBRAIC[63] =  STATES[1]*((STATES[5] -  STATES[6]*exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[65] =  CONSTANTS[52]*ALGEBRAIC[64]*(ALGEBRAIC[63]/75.0000);
ALGEBRAIC[28] = 156.500/(5.00000+exp(( - CONSTANTS[6]*ALGEBRAIC[25])/( CONSTANTS[35]*CONSTANTS[36])));
ALGEBRAIC[29] = 156.500 -  5.00000*ALGEBRAIC[28];
ALGEBRAIC[53] =  STATES[1]*((ALGEBRAIC[28] -  ALGEBRAIC[29]*exp(( CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*STATES[1])/( CONSTANTS[35]*CONSTANTS[36]))));
ALGEBRAIC[52] = 1.00000/(1.00000+exp((STATES[1] - CONSTANTS[23])/CONSTANTS[26]));
ALGEBRAIC[54] =  CONSTANTS[37]*pow(ALGEBRAIC[52], 4.00000);
ALGEBRAIC[55] =  CONSTANTS[49]*ALGEBRAIC[54]*(ALGEBRAIC[53]/45.0000);
ALGEBRAIC[70] = ALGEBRAIC[55]+ALGEBRAIC[60]+ALGEBRAIC[62]+ALGEBRAIC[65]+ALGEBRAIC[69];
}
