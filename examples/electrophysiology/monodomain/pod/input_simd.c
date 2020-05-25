#include <math.h>
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
{
STATES[0] = -79.974;
STATES[1] = -80.2;
STATES[2] = 5.9;
STATES[3] = 150.9;
STATES[4] = 5.9;
STATES[5] = 12.7;
STATES[6] = 133;
STATES[7] = 133;
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
STATES[23] = 1;
STATES[24] = 0;
STATES[25] = 0;
STATES[26] = 0;
STATES[27] = 0;
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
}
/* This function was created by opendihu at 13/08/2018 21:36:29.
 * Function is designed for 501 instances of the CellML problem. */
void computeCellMLRightHandSide(void *context, double t, double *states, double *rates, double *algebraics, double *parameters)
{
  double VOI = t;   /* current simulation time */
  double CONSTANTS[0];
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
CONSTANTS[49] = 0.1;
CONSTANTS[50] = 1.0;
CONSTANTS[51] = 0.45;
CONSTANTS[52] = 0.1;
CONSTANTS[53] = 0.1;
CONSTANTS[54] = 0.002;
CONSTANTS[55] = 1000;
CONSTANTS[56] = 0.2;
CONSTANTS[57] = 0.2;
CONSTANTS[58] = 4.5;
CONSTANTS[59] = -20;
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
CONSTANTS[98] = 300;
CONSTANTS[99] =  0.950000*CONSTANTS[65]* 3.14159265358979*pow(CONSTANTS[66], 2.00000);
CONSTANTS[100] =  0.0500000*CONSTANTS[65]* 3.14159265358979*pow(CONSTANTS[66], 2.00000);
CONSTANTS[101] =  0.0100000*CONSTANTS[99];
CONSTANTS[102] =  0.990000*CONSTANTS[99];
CONSTANTS[103] =  0.0100000*CONSTANTS[100];
CONSTANTS[104] =  0.990000*CONSTANTS[100];

  double ALGEBRAIC[35571];    /* 71 per instance * 501 instances */ 
  for(int i = 0; i < 501; i++)
  {
    rates[14028+i] = (((( ( CONSTANTS[98]*(states[9018+i]+states[9519+i]+states[10020+i]+states[10521+i]+states[11022+i]))*((states[14529+i] - states[14028+i])/CONSTANTS[101]) -  CONSTANTS[60]*((states[14028+i]/(states[14028+i]+CONSTANTS[61]))/CONSTANTS[101]))+ CONSTANTS[62]*((states[14529+i] - states[14028+i])/CONSTANTS[101]))+ - CONSTANTS[63]*((states[14028+i] - states[15030+i])/CONSTANTS[101]))+- ( ( CONSTANTS[70]*states[14028+i])*((CONSTANTS[72]+- states[16533+i])+- states[17535+i])+ - CONSTANTS[71]*states[16533+i]))+- ( ( CONSTANTS[78]*states[14028+i])*states[21543+i]+ - CONSTANTS[79]*states[19539+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[14529+i] = ((( - ( CONSTANTS[98]*(states[9018+i]+states[9519+i]+states[10020+i]+states[10521+i]+states[11022+i]))*((states[14529+i] - states[14028+i])/CONSTANTS[103])+ CONSTANTS[60]*((states[14028+i]/(states[14028+i]+CONSTANTS[61]))/CONSTANTS[103]))+ - CONSTANTS[62]*((states[14529+i] - states[14028+i])/CONSTANTS[103]))+ - CONSTANTS[64]*((states[14529+i] - states[15531+i])/CONSTANTS[103]))+- ( ( CONSTANTS[75]*states[14529+i])*(CONSTANTS[77] - states[18537+i])+ - CONSTANTS[76]*states[18537+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[15531+i] = ((( CONSTANTS[60]*((states[15030+i]/(states[15030+i]+CONSTANTS[61]))/CONSTANTS[104])+ - CONSTANTS[62]*((states[15531+i]+- states[15030+i])/CONSTANTS[104]))+ CONSTANTS[64]*((states[14529+i]+- states[15531+i])/CONSTANTS[104]))+- ( ( CONSTANTS[75]*states[15531+i])*(CONSTANTS[77]+- states[19038+i])+ - CONSTANTS[76]*states[19038+i])) -  (1000.00/1.00000)*( CONSTANTS[95]*( states[27054+i]*(0.00100000/1.00000)*states[15531+i] - CONSTANTS[97])*( states[27054+i]*(0.00100000/1.00000)*states[15531+i] - CONSTANTS[97]>0.00000 ? 1.00000 : 0.00000)*(0.00100000/1.00000)*states[27054+i]*states[15531+i] -  CONSTANTS[96]*states[27555+i]*(CONSTANTS[97] -  states[27054+i]*(0.00100000/1.00000)*states[15531+i])*(CONSTANTS[97] -  states[27054+i]*(0.00100000/1.00000)*states[15531+i]>0.00000 ? 1.00000 : 0.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[16533+i] =  ( CONSTANTS[70]*states[14028+i])*((CONSTANTS[72]+- states[16533+i])+- states[17535+i])+ - CONSTANTS[71]*states[16533+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[17034+i] =  ( CONSTANTS[70]*states[15030+i])*((CONSTANTS[72]+- states[17034+i])+- states[18036+i])+ - CONSTANTS[71]*states[17034+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[17535+i] =  ( CONSTANTS[73]*(CONSTANTS[72]+- states[16533+i]+- states[17535+i]))*states[22545+i]+ - CONSTANTS[74]*states[17535+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[18036+i] =  ( CONSTANTS[73]*(CONSTANTS[72]+- states[17034+i]+- states[18036+i]))*states[23046+i]+ - CONSTANTS[74]*states[18036+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[18537+i] =  ( CONSTANTS[75]*states[14529+i])*(CONSTANTS[77]+- states[18537+i])+ - CONSTANTS[76]*states[18537+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[19038+i] =  ( CONSTANTS[75]*states[15531+i])*(CONSTANTS[77]+- states[19038+i])+ - CONSTANTS[76]*states[19038+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[19539+i] = ( ( CONSTANTS[78]*states[14028+i])*states[21543+i]+ - CONSTANTS[79]*states[19539+i])+ - CONSTANTS[82]*((states[19539+i]+- states[20040+i])/CONSTANTS[101]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[20040+i] = ( ( CONSTANTS[78]*states[15030+i])*states[22044+i]+ - CONSTANTS[79]*states[20040+i])+ CONSTANTS[82]*((states[19539+i]+- states[20040+i])/CONSTANTS[102]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[20541+i] = ( ( CONSTANTS[80]*states[22545+i])*states[21543+i]+ - CONSTANTS[81]*states[20541+i])+ - CONSTANTS[82]*((states[20541+i]+- states[21042+i])/CONSTANTS[101]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[21042+i] = ( ( CONSTANTS[80]*states[23046+i])*states[22044+i]+ - CONSTANTS[81]*states[21042+i])+ CONSTANTS[82]*((states[20541+i]+- states[21042+i])/CONSTANTS[102]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[21543+i] = (- ( ( CONSTANTS[78]*states[14028+i])*states[21543+i]+ - CONSTANTS[79]*states[19539+i])+- ( ( CONSTANTS[80]*states[22545+i])*states[21543+i]+ - CONSTANTS[81]*states[20541+i]))+ - CONSTANTS[82]*((states[21543+i]+- states[22044+i])/CONSTANTS[101]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[22044+i] = (- ( ( CONSTANTS[78]*states[15030+i])*states[22044+i]+ - CONSTANTS[79]*states[20040+i])+- ( ( CONSTANTS[80]*states[23046+i])*states[22044+i]+ - CONSTANTS[81]*states[21042+i]))+ CONSTANTS[82]*((states[21543+i]+- states[22044+i])/CONSTANTS[102]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[22545+i] = (- ( ( CONSTANTS[73]*(CONSTANTS[72]+- states[16533+i]+- states[17535+i]))*states[22545+i]+ - CONSTANTS[74]*states[17535+i])+- ( ( CONSTANTS[80]*states[22545+i])*states[21543+i]+ - CONSTANTS[81]*states[20541+i]))+ - CONSTANTS[83]*((states[22545+i]+- states[23046+i])/CONSTANTS[101]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[23046+i] = (- ( ( CONSTANTS[73]*(CONSTANTS[72]+- states[17034+i]+- states[18036+i]))*states[23046+i]+ - CONSTANTS[74]*states[18036+i])+- ( ( CONSTANTS[80]*states[23046+i])*states[22044+i]+ - CONSTANTS[81]*states[21042+i]))+ CONSTANTS[83]*((states[22545+i]+- states[23046+i])/CONSTANTS[102]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[23547+i] = (( ( CONSTANTS[67]*states[15030+i])*states[16032+i]+ - CONSTANTS[68]*states[23547+i])+ - CONSTANTS[86]*states[23547+i])+ CONSTANTS[87]*states[25050+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[24549+i] = (((( CONSTANTS[67]*states[15030+i]*states[24048+i]+ - CONSTANTS[68]*states[24549+i])+ CONSTANTS[84]*states[16032+i])+ - CONSTANTS[85]*states[24549+i])+ ( - CONSTANTS[67]*states[15030+i])*states[24549+i])+ CONSTANTS[68]*states[25050+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[25050+i] = ((((( CONSTANTS[67]*states[15030+i]*states[24549+i]+ - CONSTANTS[68]*states[25050+i])+ CONSTANTS[86]*states[23547+i])+ - CONSTANTS[87]*states[25050+i])+ - CONSTANTS[88]*states[25050+i])+ CONSTANTS[89]*states[25551+i])+ CONSTANTS[92]*states[26052+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[25551+i] = (( CONSTANTS[88]*states[25050+i]+ - CONSTANTS[89]*states[25551+i])+ CONSTANTS[91]*states[26052+i])+ - CONSTANTS[90]*states[25551+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[26052+i] = ( - CONSTANTS[91]*states[26052+i]+ CONSTANTS[90]*states[25551+i])+ - CONSTANTS[92]*states[26052+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[26553+i] =  (0.00100000/1.00000)*( CONSTANTS[90]*states[25551+i] -  CONSTANTS[91]*states[26052+i])+ -1.00000*CONSTANTS[93]*states[26553+i]+ -1.00000*CONSTANTS[94]*((states[26553+i] - states[27054+i])/CONSTANTS[102]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[27054+i] =  CONSTANTS[94]*((states[26553+i] - states[27054+i])/CONSTANTS[104]) -  1.00000*( CONSTANTS[95]*( states[27054+i]*(0.00100000/1.00000)*states[15531+i] - CONSTANTS[97])*( states[27054+i]*(0.00100000/1.00000)*states[15531+i] - CONSTANTS[97]>0.00000 ? 1.00000 : 0.00000)*(0.00100000/1.00000)*states[27054+i]*states[15531+i] -  CONSTANTS[96]*states[27555+i]*(CONSTANTS[97] -  states[27054+i]*(0.00100000/1.00000)*states[15531+i])*(CONSTANTS[97] -  states[27054+i]*(0.00100000/1.00000)*states[15531+i]>0.00000 ? 1.00000 : 0.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[27555+i] =  1.00000*( CONSTANTS[95]*( states[27054+i]*(0.00100000/1.00000)*states[15531+i] - CONSTANTS[97])*( states[27054+i]*(0.00100000/1.00000)*states[15531+i] - CONSTANTS[97]>0.00000 ? 1.00000 : 0.00000)*(0.00100000/1.00000)*states[27054+i]*states[15531+i] -  CONSTANTS[96]*states[27555+i]*(CONSTANTS[97] -  states[27054+i]*(0.00100000/1.00000)*states[15531+i])*(CONSTANTS[97] -  states[27054+i]*(0.00100000/1.00000)*states[15531+i]>0.00000 ? 1.00000 : 0.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[6012+i] = CONSTANTS[69]+- states[16032+i]+- states[23547+i]+- states[24048+i]+- states[24549+i]+- states[25050+i]+- states[25551+i]+- states[26052+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[15030+i] = (((( - CONSTANTS[60]*((states[15030+i]/(states[15030+i]+CONSTANTS[61]))/CONSTANTS[102])+ CONSTANTS[62]*((states[15531+i]+- states[15030+i])/CONSTANTS[102]))+ CONSTANTS[63]*((states[14028+i] - states[15030+i])/CONSTANTS[102]))+- ((((((( CONSTANTS[67]*states[15030+i]*algebraics[6012+i]+ - CONSTANTS[68]*states[16032+i])+ CONSTANTS[67]*states[15030+i]*states[16032+i])+ - CONSTANTS[68]*states[23547+i])+ CONSTANTS[67]*states[15030+i]*states[24048+i])+ - CONSTANTS[68]*states[24549+i])+ CONSTANTS[67]*states[15030+i]*states[24549+i])+ - CONSTANTS[68]*states[25050+i]))+- ( ( CONSTANTS[70]*states[15030+i])*(CONSTANTS[72]+- states[17034+i]+- states[18036+i])+ - CONSTANTS[71]*states[17034+i]))+- ( ( CONSTANTS[78]*states[15030+i])*states[22044+i]+ - CONSTANTS[79]*states[20040+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[16032+i] = (((( ( CONSTANTS[67]*states[15030+i])*algebraics[6012+i]+ - CONSTANTS[68]*states[16032+i])+ ( - CONSTANTS[67]*states[15030+i])*states[16032+i])+ CONSTANTS[68]*states[23547+i])+ - CONSTANTS[84]*states[16032+i])+ CONSTANTS[85]*states[24549+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[24048+i] = (( ( - CONSTANTS[67]*states[15030+i])*states[24048+i]+ CONSTANTS[68]*states[24549+i])+ CONSTANTS[84]*algebraics[6012+i])+ - CONSTANTS[85]*states[24048+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[501+i] =  CONSTANTS[16]*((states[0+i] - CONSTANTS[21])/(1.00000 - exp(- ((states[0+i] - CONSTANTS[21])/CONSTANTS[32]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[7014+i] =  CONSTANTS[19]*exp(- ((states[0+i] - CONSTANTS[21])/CONSTANTS[34]));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[4008+i] =  algebraics[501+i]*(1.00000 - states[4008+i]) -  algebraics[7014+i]*states[4008+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[1002+i] = 1.00000/(1.00000+exp((states[0+i] - CONSTANTS[25])/CONSTANTS[28]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[7515+i] =  1000.00*exp(- ((states[0+i]+40.0000)/25.7500));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[4509+i] = (algebraics[1002+i] - states[4509+i])/algebraics[7515+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[2004+i] =  CONSTANTS[15]*((states[0+i] - CONSTANTS[20])/(1.00000 - exp(- ((states[0+i] - CONSTANTS[20])/CONSTANTS[31]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[8517+i] =  CONSTANTS[18]*exp(- ((states[0+i] - CONSTANTS[20])/CONSTANTS[33]));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[5010+i] =  algebraics[2004+i]*(1.00000 - states[5010+i]) -  algebraics[8517+i]*states[5010+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[1503+i] =  CONSTANTS[14]*exp(- ((states[0+i] - CONSTANTS[22])/CONSTANTS[29]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[8016+i] = CONSTANTS[17]/(1.00000+exp(- ((states[0+i] - CONSTANTS[22])/CONSTANTS[30])));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[5511+i] =  algebraics[1503+i]*(1.00000 - states[5511+i]) -  algebraics[8016+i]*states[5511+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[2505+i] = 1.00000/(1.00000+exp((states[0+i] - CONSTANTS[24])/CONSTANTS[27]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[9018+i] = 8571.00/(0.200000+ 5.65000*pow((states[0+i]+CONSTANTS[48])/100.000, 2.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[6012+i] = (algebraics[2505+i] - states[6012+i])/algebraics[9018+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[3006+i] =  CONSTANTS[16]*((states[501+i] - CONSTANTS[21])/(1.00000 - exp(- ((states[501+i] - CONSTANTS[21])/CONSTANTS[32]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[9519+i] =  CONSTANTS[19]*exp(- ((states[501+i] - CONSTANTS[21])/CONSTANTS[34]));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[6513+i] =  algebraics[3006+i]*(1.00000 - states[6513+i]) -  algebraics[9519+i]*states[6513+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[3507+i] = 1.00000/(1.00000+exp((states[501+i] - CONSTANTS[25])/CONSTANTS[28]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[10020+i] =  1.00000*exp(- ((states[501+i]+40.0000)/25.7500));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[7014+i] = (algebraics[3507+i] - states[7014+i])/algebraics[10020+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[4509+i] =  CONSTANTS[15]*((states[501+i] - CONSTANTS[20])/(1.00000 - exp(- ((states[501+i] - CONSTANTS[20])/CONSTANTS[31]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[11022+i] =  CONSTANTS[18]*exp(- ((states[501+i] - CONSTANTS[20])/CONSTANTS[33]));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[7515+i] =  algebraics[4509+i]*(1.00000 - states[7515+i]) -  algebraics[11022+i]*states[7515+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[4008+i] =  CONSTANTS[14]*exp(- ((states[501+i] - CONSTANTS[22])/CONSTANTS[29]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[10521+i] = CONSTANTS[17]/(1.00000+exp(- ((states[501+i] - CONSTANTS[22])/CONSTANTS[30])));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[8016+i] =  algebraics[4008+i]*(1.00000 - states[8016+i]) -  algebraics[10521+i]*states[8016+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[5010+i] = 1.00000/(1.00000+exp((states[501+i] - CONSTANTS[24])/CONSTANTS[27]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[11523+i] = 8571.00/(0.200000+ 5.65000*pow((states[501+i]+CONSTANTS[48])/100.000, 2.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[8517+i] = (algebraics[5010+i] - states[8517+i])/algebraics[11523+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[5511+i] =  0.500000*CONSTANTS[57]*exp((states[501+i] - CONSTANTS[59])/( 8.00000*CONSTANTS[58]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[12024+i] =  0.500000*CONSTANTS[57]*exp((CONSTANTS[59] - states[501+i])/( 8.00000*CONSTANTS[58]));
  }
  for(int i = 0; i < 501; i++)
  {
    rates[11523+i] =  - CONSTANTS[54]*states[11523+i]+ CONSTANTS[55]*states[9018+i]+ -4.00000*algebraics[5511+i]*states[11523+i]+ algebraics[12024+i]*states[12024+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[9018+i] =  CONSTANTS[54]*states[11523+i]+ - CONSTANTS[55]*states[9018+i]+( -4.00000*algebraics[5511+i]*states[9018+i])/CONSTANTS[56]+ CONSTANTS[56]*algebraics[12024+i]*states[9519+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[12024+i] =  4.00000*algebraics[5511+i]*states[11523+i]+ - algebraics[12024+i]*states[12024+i]+( - CONSTANTS[54]*states[12024+i])/CONSTANTS[56]+ CONSTANTS[56]*CONSTANTS[55]*states[9519+i]+ -3.00000*algebraics[5511+i]*states[12024+i]+ 2.00000*algebraics[12024+i]*states[12525+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[9519+i] = ( CONSTANTS[54]*states[12024+i])/CONSTANTS[56]+ - CONSTANTS[55]*CONSTANTS[56]*states[9519+i]+( 4.00000*algebraics[5511+i]*states[9018+i])/CONSTANTS[56]+ - CONSTANTS[56]*algebraics[12024+i]*states[9519+i]+( -3.00000*algebraics[5511+i]*states[9519+i])/CONSTANTS[56]+ 2.00000*CONSTANTS[56]*algebraics[12024+i]*states[10020+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[12525+i] =  3.00000*algebraics[5511+i]*states[12024+i]+ -2.00000*algebraics[12024+i]*states[12525+i]+( - CONSTANTS[54]*states[12525+i])/pow(CONSTANTS[56], 2.00000)+ pow(CONSTANTS[56], 2.00000)*CONSTANTS[55]*states[10020+i]+ -2.00000*algebraics[5511+i]*states[12525+i]+ 3.00000*algebraics[12024+i]*states[13026+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[10020+i] = ( 3.00000*algebraics[5511+i]*states[9519+i])/CONSTANTS[56]+ -2.00000*CONSTANTS[56]*algebraics[12024+i]*states[10020+i]+( CONSTANTS[54]*states[12525+i])/pow(CONSTANTS[56], 2.00000)+ - CONSTANTS[55]*pow(CONSTANTS[56], 2.00000)*states[10020+i]+( -2.00000*algebraics[5511+i]*states[10020+i])/CONSTANTS[56]+ 3.00000*CONSTANTS[56]*algebraics[12024+i]*states[10521+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[13026+i] =  2.00000*algebraics[5511+i]*states[12525+i]+ -3.00000*algebraics[12024+i]*states[13026+i]+( - CONSTANTS[54]*states[13026+i])/pow(CONSTANTS[56], 3.00000)+ CONSTANTS[55]*pow(CONSTANTS[56], 3.00000)*states[10521+i]+ - algebraics[5511+i]*states[13026+i]+ 4.00000*algebraics[12024+i]*states[13527+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[10521+i] = ( CONSTANTS[54]*states[13026+i])/pow(CONSTANTS[56], 3.00000)+ - CONSTANTS[55]*pow(CONSTANTS[56], 3.00000)*states[10521+i]+( 2.00000*algebraics[5511+i]*states[10020+i])/CONSTANTS[56]+ -3.00000*algebraics[12024+i]*CONSTANTS[56]*states[10521+i]+( - algebraics[5511+i]*states[10521+i])/CONSTANTS[56]+ 4.00000*CONSTANTS[56]*algebraics[12024+i]*states[11022+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[13527+i] =  algebraics[5511+i]*states[13026+i]+ -4.00000*algebraics[12024+i]*states[13527+i]+( - CONSTANTS[54]*states[13527+i])/pow(CONSTANTS[56], 4.00000)+ CONSTANTS[55]*pow(CONSTANTS[56], 4.00000)*states[11022+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[11022+i] = ( algebraics[5511+i]*states[10521+i])/CONSTANTS[56]+ -4.00000*CONSTANTS[56]*algebraics[12024+i]*states[11022+i]+( CONSTANTS[54]*states[13527+i])/pow(CONSTANTS[56], 4.00000)+ - CONSTANTS[55]*pow(CONSTANTS[56], 4.00000)*states[11022+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[15030+i] =  states[0+i]*((states[1503+i] -  states[2004+i]*exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[6513+i] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(states[2004+i]/states[1503+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[18537+i] =  states[2004+i]*exp( ( - CONSTANTS[41]*algebraics[6513+i])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[19038+i] =  CONSTANTS[40]*(pow(algebraics[18537+i], 2.00000)/(CONSTANTS[42]+pow(algebraics[18537+i], 2.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[19539+i] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(algebraics[18537+i], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*states[0+i]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[20040+i] =  algebraics[19038+i]*algebraics[19539+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[20541+i] =  algebraics[20040+i]*(algebraics[15030+i]>0.00000 ? 1.00000 : 0.00000)*(algebraics[15030+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[21042+i] =  ( CONSTANTS[38]*pow(states[4008+i], 4.00000))*states[4509+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[21543+i] =  algebraics[21042+i]*(algebraics[15030+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[23547+i] =  (1.00000/7.00000)*(exp(states[3507+i]/67.3000) - 1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[24048+i] = pow(1.00000+ 0.120000*exp( -0.100000*states[0+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*algebraics[23547+i]*exp(- ( states[0+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[24549+i] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/states[2004+i], 2.00000)*pow(1.00000+CONSTANTS[45]/states[2505+i], 3.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[25050+i] =  algebraics[24549+i]*algebraics[24048+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[2004+i] = (algebraics[20541+i]+algebraics[21543+i]+CONSTANTS[12]+ - 2.00000*algebraics[25050+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[5])+(states[1002+i] - states[2004+i])/CONSTANTS[10];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[22545+i] =  ( ( CONSTANTS[39]*pow(states[5010+i], 3.00000))*states[5511+i])*states[6012+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[22044+i] =  states[0+i]*((states[2505+i] -  states[3507+i]*exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[23046+i] =  algebraics[22545+i]*(algebraics[22044+i]/75.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[3507+i] = (algebraics[23046+i]+CONSTANTS[13]+ 3.00000*algebraics[25050+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[5])+(states[3006+i] - states[3507+i])/CONSTANTS[11];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[0+i] =  (1000.00/1.00000)*((states[0+i] - states[501+i])/CONSTANTS[2]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[13026+i] = 156.500/(5.00000+exp(( - CONSTANTS[6]*algebraics[6513+i])/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[13527+i] = 156.500 -  5.00000*algebraics[13026+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[17034+i] =  states[0+i]*((algebraics[13026+i] -  algebraics[13527+i]*exp(( CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[16533+i] = 1.00000/(1.00000+exp((states[0+i] - CONSTANTS[23])/CONSTANTS[26]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[17535+i] =  CONSTANTS[37]*pow(algebraics[16533+i], 4.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[18036+i] =  algebraics[17535+i]*(algebraics[17034+i]/45.0000);
  }
  /* explicit parameter */
  /* ALGEBRAIC[32] = (VOI>=0.00000&&VOI<0.500000 ? 150.000 : VOI>=50.0000&&VOI<50.5000 ? 150.000 : VOI>=100.000&&VOI<100.500 ? 150.000 : VOI>=150.000&&VOI<150.500 ? 150.000 : VOI>=200.000&&VOI<200.500 ? 150.000 : VOI>=250.000&&VOI<250.500 ? 150.000 : VOI>=300.000&&VOI<300.500 ? 150.000 : VOI>=350.000&&VOI<350.500 ? 150.000 : VOI>=400.000&&VOI<400.500 ? 150.000 : 0.00000);*/
  for(int i = 0; i < 501; i++)
  {
    algebraics[25551+i] = algebraics[18036+i]+algebraics[20541+i]+algebraics[21543+i]+algebraics[23046+i]+algebraics[25050+i]+- parameters[0+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[0+i] = - ((algebraics[25551+i]+algebraics[0+i])/CONSTANTS[0]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[15531+i] =  states[501+i]*((states[1503+i] -  states[1002+i]*exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[12525+i] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(states[1002+i]/states[1503+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[28056+i] =  states[1002+i]*exp( ( - CONSTANTS[41]*algebraics[12525+i])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[28557+i] =  CONSTANTS[40]*(pow(algebraics[28056+i], 2.00000)/(CONSTANTS[42]+pow(algebraics[28056+i], 2.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[29058+i] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(algebraics[28056+i], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*states[501+i]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[29559+i] =  algebraics[28557+i]*algebraics[29058+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[30060+i] =  CONSTANTS[50]*algebraics[29559+i]*(algebraics[15531+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[30561+i] =  ( CONSTANTS[38]*pow(states[6513+i], 4.00000))*states[7014+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[31062+i] =  CONSTANTS[51]*algebraics[30561+i]*(algebraics[15531+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[33066+i] =  (1.00000/7.00000)*(exp(states[3006+i]/67.3000) - 1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[33567+i] = pow(1.00000+ 0.120000*exp( -0.100000*states[501+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*algebraics[33066+i]*exp(- ( states[501+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[34068+i] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/states[1002+i], 2.00000)*pow(1.00000+CONSTANTS[45]/states[2505+i], 3.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[34569+i] =  CONSTANTS[53]*algebraics[34068+i]*algebraics[33567+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[1503+i] =  - CONSTANTS[9]*((algebraics[30060+i]+algebraics[31062+i]+CONSTANTS[12]+ - 2.00000*algebraics[34569+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3])) - (algebraics[20541+i]+algebraics[21543+i]+CONSTANTS[12]+ -2.00000*algebraics[25050+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[4]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[1002+i] = (algebraics[30060+i]+algebraics[31062+i]+CONSTANTS[12]+ - 2.00000*algebraics[34569+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3]) - (states[1002+i] - states[2004+i])/CONSTANTS[7];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[32064+i] =  ( ( CONSTANTS[39]*pow(states[7515+i], 3.00000))*states[8016+i])*states[8517+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[31563+i] =  states[501+i]*((states[2505+i] -  states[3006+i]*exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[32565+i] =  CONSTANTS[52]*algebraics[32064+i]*(algebraics[31563+i]/75.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[2505+i] =  - CONSTANTS[9]*((algebraics[32565+i]+CONSTANTS[13]+ 3.00000*algebraics[34569+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3])) - (algebraics[23046+i]+CONSTANTS[13]+ 3.00000*algebraics[25050+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[4]);
  }
  for(int i = 0; i < 501; i++)
  {
    rates[3006+i] = (algebraics[32565+i]+CONSTANTS[13]+ 3.00000*algebraics[34569+i])/( (1000.00/1.00000)*CONSTANTS[6]*CONSTANTS[3]) - (states[3006+i] - states[3507+i])/CONSTANTS[8];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[14028+i] = 156.500/(5.00000+exp(( - CONSTANTS[6]*algebraics[12525+i])/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[14529+i] = 156.500 -  5.00000*algebraics[14028+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[26553+i] =  states[501+i]*((algebraics[14028+i] -  algebraics[14529+i]*exp(( CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[26052+i] = 1.00000/(1.00000+exp((states[501+i] - CONSTANTS[23])/CONSTANTS[26]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[27054+i] =  CONSTANTS[37]*pow(algebraics[26052+i], 4.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[27555+i] =  CONSTANTS[49]*algebraics[27054+i]*(algebraics[26553+i]/45.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[35070+i] = algebraics[27555+i]+algebraics[30060+i]+algebraics[31062+i]+algebraics[32565+i]+algebraics[34569+i];
  }
  for(int i = 0; i < 501; i++)
  {
    rates[501+i] = - ((algebraics[35070+i] - algebraics[0+i]/CONSTANTS[1])/CONSTANTS[0]);
  }
}
void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
  for(int i = 0; i < 501; i++)
  {
    algebraics[6012+i] = CONSTANTS[69]+- states[16032+i]+- states[23547+i]+- states[24048+i]+- states[24549+i]+- states[25050+i]+- states[25551+i]+- states[26052+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[501+i] =  CONSTANTS[16]*((states[0+i] - CONSTANTS[21])/(1.00000 - exp(- ((states[0+i] - CONSTANTS[21])/CONSTANTS[32]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[7014+i] =  CONSTANTS[19]*exp(- ((states[0+i] - CONSTANTS[21])/CONSTANTS[34]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[1002+i] = 1.00000/(1.00000+exp((states[0+i] - CONSTANTS[25])/CONSTANTS[28]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[7515+i] =  1000.00*exp(- ((states[0+i]+40.0000)/25.7500));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[2004+i] =  CONSTANTS[15]*((states[0+i] - CONSTANTS[20])/(1.00000 - exp(- ((states[0+i] - CONSTANTS[20])/CONSTANTS[31]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[8517+i] =  CONSTANTS[18]*exp(- ((states[0+i] - CONSTANTS[20])/CONSTANTS[33]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[1503+i] =  CONSTANTS[14]*exp(- ((states[0+i] - CONSTANTS[22])/CONSTANTS[29]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[8016+i] = CONSTANTS[17]/(1.00000+exp(- ((states[0+i] - CONSTANTS[22])/CONSTANTS[30])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[2505+i] = 1.00000/(1.00000+exp((states[0+i] - CONSTANTS[24])/CONSTANTS[27]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[9018+i] = 8571.00/(0.200000+ 5.65000*pow((states[0+i]+CONSTANTS[48])/100.000, 2.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[3006+i] =  CONSTANTS[16]*((states[501+i] - CONSTANTS[21])/(1.00000 - exp(- ((states[501+i] - CONSTANTS[21])/CONSTANTS[32]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[9519+i] =  CONSTANTS[19]*exp(- ((states[501+i] - CONSTANTS[21])/CONSTANTS[34]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[3507+i] = 1.00000/(1.00000+exp((states[501+i] - CONSTANTS[25])/CONSTANTS[28]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[10020+i] =  1.00000*exp(- ((states[501+i]+40.0000)/25.7500));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[4509+i] =  CONSTANTS[15]*((states[501+i] - CONSTANTS[20])/(1.00000 - exp(- ((states[501+i] - CONSTANTS[20])/CONSTANTS[31]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[11022+i] =  CONSTANTS[18]*exp(- ((states[501+i] - CONSTANTS[20])/CONSTANTS[33]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[4008+i] =  CONSTANTS[14]*exp(- ((states[501+i] - CONSTANTS[22])/CONSTANTS[29]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[10521+i] = CONSTANTS[17]/(1.00000+exp(- ((states[501+i] - CONSTANTS[22])/CONSTANTS[30])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[5010+i] = 1.00000/(1.00000+exp((states[501+i] - CONSTANTS[24])/CONSTANTS[27]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[11523+i] = 8571.00/(0.200000+ 5.65000*pow((states[501+i]+CONSTANTS[48])/100.000, 2.00000));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[5511+i] =  0.500000*CONSTANTS[57]*exp((states[501+i] - CONSTANTS[59])/( 8.00000*CONSTANTS[58]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[12024+i] =  0.500000*CONSTANTS[57]*exp((CONSTANTS[59] - states[501+i])/( 8.00000*CONSTANTS[58]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[15030+i] =  states[0+i]*((states[1503+i] -  states[2004+i]*exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[6513+i] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(states[2004+i]/states[1503+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[18537+i] =  states[2004+i]*exp( ( - CONSTANTS[41]*algebraics[6513+i])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[19038+i] =  CONSTANTS[40]*(pow(algebraics[18537+i], 2.00000)/(CONSTANTS[42]+pow(algebraics[18537+i], 2.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[19539+i] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(algebraics[18537+i], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*states[0+i]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[20040+i] =  algebraics[19038+i]*algebraics[19539+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[20541+i] =  algebraics[20040+i]*(algebraics[15030+i]>0.00000 ? 1.00000 : 0.00000)*(algebraics[15030+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[21042+i] =  ( CONSTANTS[38]*pow(states[4008+i], 4.00000))*states[4509+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[21543+i] =  algebraics[21042+i]*(algebraics[15030+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[23547+i] =  (1.00000/7.00000)*(exp(states[3507+i]/67.3000) - 1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[24048+i] = pow(1.00000+ 0.120000*exp( -0.100000*states[0+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*algebraics[23547+i]*exp(- ( states[0+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[24549+i] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/states[2004+i], 2.00000)*pow(1.00000+CONSTANTS[45]/states[2505+i], 3.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[25050+i] =  algebraics[24549+i]*algebraics[24048+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[22545+i] =  ( ( CONSTANTS[39]*pow(states[5010+i], 3.00000))*states[5511+i])*states[6012+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[22044+i] =  states[0+i]*((states[2505+i] -  states[3507+i]*exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[23046+i] =  algebraics[22545+i]*(algebraics[22044+i]/75.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[0+i] =  (1000.00/1.00000)*((states[0+i] - states[501+i])/CONSTANTS[2]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[13026+i] = 156.500/(5.00000+exp(( - CONSTANTS[6]*algebraics[6513+i])/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[13527+i] = 156.500 -  5.00000*algebraics[13026+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[17034+i] =  states[0+i]*((algebraics[13026+i] -  algebraics[13527+i]*exp(( CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*states[0+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[16533+i] = 1.00000/(1.00000+exp((states[0+i] - CONSTANTS[23])/CONSTANTS[26]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[17535+i] =  CONSTANTS[37]*pow(algebraics[16533+i], 4.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[18036+i] =  algebraics[17535+i]*(algebraics[17034+i]/45.0000);
  }
  /* explicit parameter */
  /* ALGEBRAIC[32] = (VOI>=0.00000&&VOI<0.500000 ? 150.000 : VOI>=50.0000&&VOI<50.5000 ? 150.000 : VOI>=100.000&&VOI<100.500 ? 150.000 : VOI>=150.000&&VOI<150.500 ? 150.000 : VOI>=200.000&&VOI<200.500 ? 150.000 : VOI>=250.000&&VOI<250.500 ? 150.000 : VOI>=300.000&&VOI<300.500 ? 150.000 : VOI>=350.000&&VOI<350.500 ? 150.000 : VOI>=400.000&&VOI<400.500 ? 150.000 : 0.00000);*/
  for(int i = 0; i < 501; i++)
  {
    algebraics[25551+i] = algebraics[18036+i]+algebraics[20541+i]+algebraics[21543+i]+algebraics[23046+i]+algebraics[25050+i]+- parameters[0+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[15531+i] =  states[501+i]*((states[1503+i] -  states[1002+i]*exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[12525+i] =  (( CONSTANTS[35]*CONSTANTS[36])/CONSTANTS[6])*log(states[1002+i]/states[1503+i]);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[28056+i] =  states[1002+i]*exp( ( - CONSTANTS[41]*algebraics[12525+i])*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[28557+i] =  CONSTANTS[40]*(pow(algebraics[28056+i], 2.00000)/(CONSTANTS[42]+pow(algebraics[28056+i], 2.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[29058+i] = 1.00000 - pow(1.00000+( CONSTANTS[43]*(1.00000+pow(algebraics[28056+i], 2.00000)/CONSTANTS[42]))/( pow(CONSTANTS[46], 2.00000)*exp(( 2.00000*(1.00000 - CONSTANTS[41])*states[501+i]*CONSTANTS[6])/( CONSTANTS[35]*CONSTANTS[36]))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[29559+i] =  algebraics[28557+i]*algebraics[29058+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[30060+i] =  CONSTANTS[50]*algebraics[29559+i]*(algebraics[15531+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[30561+i] =  ( CONSTANTS[38]*pow(states[6513+i], 4.00000))*states[7014+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[31062+i] =  CONSTANTS[51]*algebraics[30561+i]*(algebraics[15531+i]/50.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[33066+i] =  (1.00000/7.00000)*(exp(states[3006+i]/67.3000) - 1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[33567+i] = pow(1.00000+ 0.120000*exp( -0.100000*states[501+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))+ 0.0400000*algebraics[33066+i]*exp(- ( states[501+i]*(CONSTANTS[6]/( CONSTANTS[35]*CONSTANTS[36])))), -1.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[34068+i] =  CONSTANTS[6]*(CONSTANTS[47]/( pow(1.00000+CONSTANTS[44]/states[1002+i], 2.00000)*pow(1.00000+CONSTANTS[45]/states[2505+i], 3.00000)));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[34569+i] =  CONSTANTS[53]*algebraics[34068+i]*algebraics[33567+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[32064+i] =  ( ( CONSTANTS[39]*pow(states[7515+i], 3.00000))*states[8016+i])*states[8517+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[31563+i] =  states[501+i]*((states[2505+i] -  states[3006+i]*exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( -1.00000*CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[32565+i] =  CONSTANTS[52]*algebraics[32064+i]*(algebraics[31563+i]/75.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[14028+i] = 156.500/(5.00000+exp(( - CONSTANTS[6]*algebraics[12525+i])/( CONSTANTS[35]*CONSTANTS[36])));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[14529+i] = 156.500 -  5.00000*algebraics[14028+i];
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[26553+i] =  states[501+i]*((algebraics[14028+i] -  algebraics[14529+i]*exp(( CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36])))/(1.00000 - exp(( CONSTANTS[6]*states[501+i])/( CONSTANTS[35]*CONSTANTS[36]))));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[26052+i] = 1.00000/(1.00000+exp((states[501+i] - CONSTANTS[23])/CONSTANTS[26]));
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[27054+i] =  CONSTANTS[37]*pow(algebraics[26052+i], 4.00000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[27555+i] =  CONSTANTS[49]*algebraics[27054+i]*(algebraics[26553+i]/45.0000);
  }
  for(int i = 0; i < 501; i++)
  {
    algebraics[35070+i] = algebraics[27555+i]+algebraics[30060+i]+algebraics[31062+i]+algebraics[32565+i]+algebraics[34569+i];
  }
}

