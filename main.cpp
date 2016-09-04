
#include"header.h"

using namespace std;

int main()
{
  ifstream QHfile;

  double step_size = 1e-4;
  double bobyqa_stop_criteria = 1e-4;
// step_size = std::sqrt( std::sqrt(epsilonGet<double>()) );   // step_size = 0.0001;
  p_step_size = &step_size;
 //  Out<<setprecision(6)<<endl;
  // cout<<"step_size auto calculated = "<<step_size<<endl;

    int nsub;
    int mQHpts, i, j;
    mQHpts = 1; dbmat aQHPts1(mQHpts,2);
    aQHPts1 = 0 , 1.77245385090552;

    mQHpts = 3; dbmat aQHPts3(mQHpts,2);
    aQHPts3 = -1.22474487139159, 0.29540897515092,
-1.46354448532132e-16, 1.18163590060368,
1.22474487139159, 0.29540897515092;

    mQHpts = 5; dbmat aQHPts5(mQHpts,2); aQHPts5 = -2.02018287045609,0.0199532420590459,
-0.958572464613819, 0.393619323152241,
2.73349801485409e-17, 0.945308720482941,
0.958572464613819, 0.393619323152241,
2.02018287045609, 0.0199532420590459;

    mQHpts = 7; dbmat aQHPts7(mQHpts,2);
    aQHPts7 = -2.65196135683523, 0.000971781245099518,
-1.67355162876747, 0.0545155828191271,
-0.816287882858965, 0.425607252610128,
8.15487421110398e-17, 0.810264617556807,
0.816287882858964, 0.425607252610128,
1.67355162876747, 0.0545155828191269,
2.65196135683523, 0.000971781245099514;

    mQHpts = 9; dbmat aQHPts9(mQHpts,2);
    aQHPts9 = -3.19099320178153, 3.96069772632644e-05,
-2.26658058453184, 0.00494362427553695,
-1.46855328921667, 0.0884745273943766,
-0.723551018752838, 0.432651559002556,
-4.02438494023542e-17, 0.720235215606051,
0.723551018752837, 0.432651559002556,
1.46855328921667, 0.0884745273943768,
2.26658058453184, 0.00494362427553694,
3.19099320178153, 3.96069772632643e-05;

    mQHpts = 11; dbmat aQHPts11(mQHpts,2);
    aQHPts11 = -3.66847084655958, 1.43956039371426e-06,
-2.78329009978165, 0.000346819466323345,
-2.02594801582576, 0.0119113954449115,
-1.32655708449493, 0.117227875167708,
-0.6568095668821, 0.429359752356125,
2.50185338164225e-16, 0.654759286914591,
0.6568095668821, 0.429359752356125,
1.32655708449493, 0.117227875167708,
2.02594801582576, 0.0119113954449115,
2.78329009978165, 0.000346819466323345,
3.66847084655958, 1.43956039371426e-06;

    mQHpts = 13; dbmat aQHPts13(mQHpts,2);
    aQHPts13 =  -4.10133759617864, 4.82573185007313e-08,
-3.24660897837241, 2.04303604027071e-05,
-2.51973568567824, 0.00120745999271938,
-1.85310765160151, 0.0208627752961699,
-1.22005503659075, 0.140323320687023,
-0.60576387917106, 0.421616296898543,
1.48984508284312e-16, 0.604393187921162,
0.60576387917106, 0.421616296898543,
1.22005503659075, 0.140323320687023,
1.85310765160151, 0.02086277529617,
2.51973568567824, 0.00120745999271939,
3.24660897837241,2.04303604027071e-05,
4.10133759617864, 4.82573185007313e-08;

  /*
    mQHpts = 15; dbmat aQHPts15(mQHpts,2);
    aQHPts15 =  -4.49999070730939, 1.52247580425352e-09,
-3.66995037340445 1.05911554771106e-06,
-2.9671669279056 0.0001000044412325,
-2.32573248617386 0.00277806884291276,
-1.71999257518649 0.0307800338725461,
-1.13611558521092 0.158488915795936,
-0.565069583255576 0.412028687498898,
-3.84143836181876e-16 0.564100308726418,
0.565069583255576 0.412028687498898,
1.13611558521092 0.158488915795935
1.71999257518649 0.030780033872546
2.32573248617386 0.00277806884291278
2.9671669279056 0.0001000044412325
3.66995037340445 1.05911554771107e-06
4.49999070730939 1.52247580425352e-09;

     mQHpts = 17; dbmat aQHPts15(mQHpts,2);
    aQHPts17 =  -4.87134519367441 4.58057893079865e-11
-4.06194667587548 4.97707898163074e-08
-3.37893209114149 7.11228914002128e-06
-2.75776291570389 0.000298643286697752
-2.17350282666662 0.00506734995762754
-1.61292431422123 0.0409200341497563
-1.06764872574345 0.172648297670097
-0.531633001342655 0.401826469470412
3.92660474520814e-16 0.530917937624865
0.531633001342655 0.401826469470412
1.06764872574345 0.172648297670098
1.61292431422123 0.0409200341497562
2.17350282666662 0.00506734995762753
2.75776291570389 0.000298643286697751
3.37893209114149 7.11228914002132e-06
4.06194667587548 4.97707898163076e-08
4.8713451936744 4.58057893079866e-11;

    mQHpts = 19; dbmat aQHPts19(mQHpts,2);
    aQHPts19 = 5.22027169053748, 1.32629709449851e-12,
-4.42853280660378, 2.16305100986355e-09,
-3.76218735196402, 4.48824314722314e-07,
-3.1578488183476, 2.72091977631616e-05,
-2.59113378979454, 0.000670877521407184,
-2.04923170985062, 0.00798886677772299,
-1.52417061939353, 0.0508103869090521,
-1.01036838713431, 0.183632701306997,
-0.503520163423888, 0.39160898861303,
-9.75022669146959e-17, 0.502974888276186,
0.503520163423888, 0.39160898861303,
1.01036838713431, 0.183632701306997,
1.52417061939353, 0.0508103869090522,
2.04923170985062, 0.00798886677772298,
2.59113378979454, 0.000670877521407181,
3.1578488183476, 2.72091977631617e-05,
3.76218735196402, 4.48824314722307e-07,
4.42853280660378, 2.16305100986356e-09,
5.22027169053748, 1.32629709449851e-12;
   */

    ifstream in1;
    in1.open("simulation_configurations.txt");
    if(!in1.is_open()) { std::cout<<"fail to open simulation configurations file\n"; return 0; }
    string str; string seed0;
    string main_out_file_name, log_file_name, dat_file_name;
    int num_simu_replication;
    dbcolvec fix_ef_glmm_true(4), fix_ef_lme_true(4);
    dbcolvec var_components(4);
    int ni_per_sub;


    for(i = 1; i < 30; ++i)
    {
        switch (i)
        {
        case 2: //  integers
            getline(in1,str);
            num_simu_replication = atoi(str.c_str());
            break;

        case 5: //  integers
            getline(in1,str);
            nsub = atoi(str.c_str());
            break;

        case 8: // string
            getline(in1,seed0);
            break;

          case 11:  // vectors
            for(j=0; j < (4-1); ++j)
            {
                getline(in1,str,',');
                fix_ef_glmm_true(j) = strtod(str.c_str(),NULL);
            }
            getline(in1,str);
            fix_ef_glmm_true(4-1) = strtod(str.c_str(),NULL);
            break;

         case 14:  // vectors
          for(j=0; j < (4-1); ++j)
            {
                getline(in1,str,',');
                fix_ef_lme_true(j) = strtod(str.c_str(),NULL);
 // std::cout<<"fix_ef_lme_true(j)\n"<<fix_ef_lme_true(j)<<"\n";
            }
            getline(in1,str);
            fix_ef_lme_true(4-1) = strtod(str.c_str(),NULL);
            break;

         case 17:  // vectors
          for(j=0; j < (4-1); ++j)
            {
                getline(in1,str,',');
                var_components(j) = strtod(str.c_str(),NULL);
            }
            getline(in1,str);
            var_components(4-1) = strtod(str.c_str(),NULL);
            break;

     case 20: //  integers
            getline(in1,str);
            ni_per_sub = atoi(str.c_str());
            break;

    case 23: // string
     getline(in1,main_out_file_name);
     break;

    case 26: // string
     getline(in1,dat_file_name);
     break;

       case 29: // string
     getline(in1,log_file_name);
     break;

         default:
            getline(in1,str);
        }
    }


  ofstream Out; Out.open(main_out_file_name.c_str());
  ofstream out_dat; out_dat.open(dat_file_name.c_str());
  ofstream out_log; out_log.open(log_file_name.c_str());


   Out<<"seed="<<seed0<<"\n";
   Out<<"nsub="<<(1*nsub)<<"\n";
  Out<<"ni_per_sub\n"<<1*ni_per_sub<<"\n";
   Out<<"num_simu_replication="<<(1*num_simu_replication)<<"\n";
   Out<<"var_components: var(bi_log_odds), var(bi_amount_vomit),corr(bi_log_odds,bi_amount_vomit), var(e_ij)\n"<<1*trans(var_components);
    Out<<"fix_ef_glmm_true, 3rd for HTi, 4rd for LTi \n"<<1*trans(fix_ef_glmm_true);
    Out<<"fix_ef_lme_true, 3rd for HTi, 4rd for LTi \n"<<1*trans(fix_ef_lme_true);

   int num_chi_bar_sam = 10000;
   dbmat L_cov_mle_chi_bar_cdf_simu(4,4);
   dbcolvec proj_of_mle_chi_bar_cdf_simu(4), std_mvn_sam(4), L_mvn_sam(4);
   dbcolvec lb_chi_bar_cdf_simu(4), ub_chi_bar_cdf_simu(4);
   lb_chi_bar_cdf_simu = 0;
   ub_chi_bar_cdf_simu = BIGNUM;
   dbcolvec chi_bar_weights(5);


  int ni = ni_per_sub, nsub3 = 3*nsub, total_ni = nsub3*ni;
  int ncol = 6; // dat: id, u_ij, v_ij, t_ij, HT_i, LT_i, first 1/3 Control group,
  dbmat dat(total_ni,ncol);

  dbcolvec ni_vec(nsub3);
  ni_vec = ni_per_sub;

  // Out<<"ni_vec\n"<<trans(ni_vec);

   // para_vec: 0-3 fixed effects for log-odds-omit; 4-7 fixed effects for amount-omit; 8:log_sd_eij; 9-11 eta determining cov(bi0,bi1)

   int num_para = 12, dim_ranf = 2, dim_eta = dim_ranf*(dim_ranf+1)/2;
   dbcolvec para_vec_true(num_para), para_vec_true_ord(num_para);
   dbmat cov_ranef_true(dim_ranf,dim_ranf), cov_ranef_true_inv(dim_ranf,dim_ranf);
   dbmat cov_ranef_L_true(dim_ranf,dim_ranf), cov_ranef_inv_L(dim_ranf,dim_ranf);
   dbcolvec cov_ranef_Lvec_true(dim_eta), cov_ranef_inv_Lvec_true(dim_eta);
   set_subm(para_vec_true,range(0,3),range(0,3)) = fix_ef_glmm_true;
   set_subm(para_vec_true,range(4,7),range(0,3)) = fix_ef_lme_true;

   para_vec_true(8) = std::log(std::sqrt(var_components(3)));
   cov_ranef_true(0,0) = var_components(0);
   cov_ranef_true(1,1) = var_components(1);
   cov_ranef_true(0,1) = var_components(2)*std::sqrt(cov_ranef_true(0,0)* cov_ranef_true(1,1));
   cov_ranef_true(1,0) = cov_ranef_true(0,1);
   // Out<<"cov_ranef_true\n"<<cov_ranef_true;

    cov_ranef_true_inv = inv(cov_ranef_true);

  // Out<<"2*cov_ranef_true="<<2*cov_ranef_true<<"\n";
 //   Out<<"cov_ranef_true="<<cov_ranef_true<<"\n";
  //  Out<<"det(cov_ranef_true) ="<<det(cov_ranef_true)<<"\n";
   cov_ranef_L_true = chol(cov_ranef_true);
   cov_ranef_inv_L = chol(cov_ranef_true_inv);
  //  Out<<"cov_ranef_L_true ="<<cov_ranef_L_true<<"\n";

   LToEta(cov_ranef_L_true,dim_ranf,cov_ranef_Lvec_true);
   LToEta(cov_ranef_inv_L,dim_ranf,cov_ranef_inv_Lvec_true);

  // Out<<"cov_ranef_inv_L ="<<trans(cov_ranef_inv_L)<<"\n";
   set_subm(para_vec_true,range(9,11),range(0,0)) = cov_ranef_Lvec_true;

  //Out<<"para_vec_true\n"<<trans(para_vec_true);

    dbcolvec para_vec_true_Lap(num_para);
    para_vec_true_Lap = para_vec_true;
    set_subm(para_vec_true_Lap,range(9,11),range(0,0)) = cov_ranef_inv_Lvec_true;

   //*******************************************************
   // ord parameterization is such that theta_ord <=0 reflects order constraints
   //**** ***************************************************
    para_vec_true_ord = para_vec_true;
   para_vec_true_ord(2) = para_vec_true(2) - para_vec_true(3);
   para_vec_true_ord(6) = para_vec_true(6) - para_vec_true(7);
   //Out<<"para_vec_true_ord\n"<<trans(para_vec_true_ord);

   //para_vec_simu_ord is HT-LT parameterization
   dbcolvec para_vec_simu(num_para),para_vec_simu_ord(num_para) , para_vec_simu_Lap(num_para);
   dlib::rand rd; rd.set_seed(seed0);
   int Two_num_para_plus1 = 2*num_para+1;
   dbcolvec lb(num_para), ub(num_para), ub_ord(num_para);
   lb = -BIGNUM; ub = BIGNUM; ub_ord = BIGNUM;
   ub_ord(2) = 0; ub_ord(3) = 0;
   ub_ord(6) = 0; ub_ord(7) = 0;
   double min_neg_loglike_unc;
   dbcolvec grad_at_unc_mle(num_para), infLoglike_para_vec(num_para);
   dbmat jm_hess(num_para,num_para), cov_jm_mle_unc(num_para,num_para);
   dbmat cov_mle_HT_LT(4,4), inv_cov_mle_HT_LT(4,4);
   dbcolvec unc_mle_HT_LT(4), con_mle_HT_LT(4);
   double jm_T_wald_ts, jm_T_wald_os, pva_jm_T_wald_ts, pva_jm_T_wald_os;
   int prop_rec_jm_T_wald_ts = 0, prop_rec_jm_T_wald_os = 0;

   dbcolvec se_unc_mle_HT_LT(4);
   dbcolvec t_ratio_unc_mle_HT_LT(4),t_ratio_con_mle_HT_LT(4);
   dbcolvec uni_p_val_unc_mle_HT_LT(4),uni_p_val_con_mle_HT_LT(4);
   int s1_para;


   dbcolvec para_vec_glmm(5), para_vec_glmm_true(5), para_vec_glmm_lb(5), para_vec_glmm_ub(5);
  para_vec_glmm_lb = -BIGNUM; para_vec_glmm_ub = BIGNUM;
  set_subm(para_vec_glmm_true,range(0,3),range(0,0)) = subm(para_vec_true,range(0,3),range(0,0));
   para_vec_glmm_true(4) = std::log(cov_ranef_true(0,0));
   dbmat glmm_hess(5,5), cov_glmm_mle_unc(5,5);

  dbcolvec para_vec_lme(6), para_vec_lme_true(6), para_vec_lme_lb(6), para_vec_lme_ub(6);
  para_vec_lme_lb = -BIGNUM; para_vec_lme_ub = BIGNUM;
  set_subm(para_vec_lme_true,range(0,3),range(0,0)) = subm(para_vec_true,range(4,7),range(0,0));
   para_vec_lme_true(4) = std::log(cov_ranef_true(1,1));
   para_vec_lme_true(5) = para_vec_true(8);
   dbmat lme_hess(6,6), cov_lme_mle_unc(6,6);


   dbcolvec para_vec_glmm_ord(5), para_vec_glmm_ord_true(5), para_vec_glmm_ord_lb(5), para_vec_glmm_ord_ub(5);
   para_vec_glmm_ord_lb = -BIGNUM; para_vec_glmm_ord_ub = BIGNUM;
   para_vec_glmm_ord_ub(2) = 0;    para_vec_glmm_ord_ub(3) = 0;
   para_vec_glmm_ord_true = para_vec_glmm_true;
   para_vec_glmm_ord_true(2) = para_vec_glmm_true(2) - para_vec_glmm_true(3); // HT- LT
   para_vec_glmm_ord_true(4) = std::log(cov_ranef_true(0,0));


   dbcolvec para_vec_lme_ord(6), para_vec_lme_ord_true(6), para_vec_lme_ord_lb(6), para_vec_lme_ord_ub(6);
   para_vec_lme_ord_lb = -BIGNUM; para_vec_lme_ord_ub = BIGNUM;
   para_vec_lme_ord_ub(2) = 0;    para_vec_lme_ord_ub(3) = 0;
      para_vec_lme_ord_true = para_vec_lme_true;
     para_vec_lme_ord_true(2) = para_vec_lme_true(2) - para_vec_lme_true(3);
   para_vec_lme_ord_true(4) = std::log(cov_ranef_true(1,1));
   para_vec_lme_ord_true(5) = para_vec_true(8);


    dbcolvec para_vec_lme_true_Lap(6), para_vec_lme_Lap(6);
    para_vec_lme_true_Lap = para_vec_lme_true;
    set_subm(para_vec_lme_true_Lap,range(4,4),range(0,0)) = std::log(std::sqrt(1/cov_ranef_true(1,1)));


   p_GH_points = &aQHPts9;
   int ConvFail;
   dbcolvec unc_mle_beta_o(2), con_mle_beta_o(2), unc_mle_beta_v(2), con_mle_beta_v(2);
   dbmat cov_beta_o(2,2), cov_beta_v(2,2);
   double glmm_T_wald_ts, lme_T_wald_ts, glmm_T_wald_os, lme_T_wald_os;
   double pva_glmm_T_wald_ts, pva_lme_T_wald_ts, pva_glmm_T_wald_os, pva_lme_T_wald_os;
   double t_ratio_unc_beta_o_HTi_min_LTi, t_ratio_unc_beta_o_LTi, t_ratio_unc_beta_v_HTi_min_LTi, t_ratio_unc_beta_v_LTi;
   int prop_rec_glmm_T_wald_ts = 0, prop_rec_glmm_T_wald_os = 0, prop_rec_glmm_T_lr_ts = 0, prop_rec_glmm_T_lr_os = 0;
   int prop_rec_lme_T_wald_ts = 0, prop_rec_lme_T_wald_os = 0, prop_rec_lme_T_lr_ts = 0, prop_rec_lme_T_lr_os = 0;
   // out_log<<"con_beta_o_HTi_min_LTi"<<","<<"unc_beta_o_HTi_min_LTi"<<","<<"t_ratio_unc_beta_o_HTi_min_LTi"<<","<<"con_beta_o_LTi"<<","<<"unc_beta_o_LTi"<<","<<"t_ratio_unc_beta_o_LTi"<<","<<"glmm_T_wald_ts"<<","<<"glmm_T_wald_os"<<","<<"pva_glmm_T_wald_ts"<<","<<"pva_glmm_T_wald_os"<<","<<"con_beta_v_HTi_min_LTi"<<","<<"unc_beta_v_HTi_min_LTi"<<","<<"t_ratio_unc_beta_v_HTi_min_LTi"<<","<<"con_beta_v_LTi"<<","<<"unc_beta_v_LTi"<<","<<"t_ratio_unc_beta_v_LTi"<<","<<"lme_T_wald_ts"<<","<<"lme_T_wald_os"<<","<<"pva_lme_T_wald_ts"<<"pva_lme_T_wald_os";

     out_log<<"unc_jm_est_HT_v"<<","<<"unc_jm_est_LT_v"<<","<<"unc_jm_est_HT_o"<<","<<"unc_jm_est_LT_o"<<","<<"con_jm_est_HT_v"<<","<<"con_jm_est_LT_v"<<","<<"con_jm_est_HT_o"<<","<<"con_jm_est_LT_o"<<","<<"pva_jm_T_wald_os"<<","<<"pva_jm_T_wald_ts"<<","<<"uni_p_val_unc_mle_HT_LT_0"<<","<<"uni_p_val_unc_mle_HT_LT_1"<<","<<"uni_p_val_unc_mle_HT_LT_2"<<","<<"uni_p_val_unc_mle_HT_LT_3"<<","<<"uni_p_val_con_mle_HT_LT_0"<<","<<"uni_p_val_con_mle_HT_LT_1"<<","<<"uni_p_val_con_mle_HT_LT_2"<<","<<"uni_p_val_con_mle_HT_LT_3"<<"\n";

   double unc_jm_est_HT_o, mean_unc_jm_est_HT_o = 0, mean_squared_unc_jm_est_HT_o = 0, true_HT_o, con_jm_est_HT_o, mean_con_jm_est_HT_o = 0, mean_squared_con_jm_est_HT_o = 0;

   double unc_jm_est_LT_o, mean_unc_jm_est_LT_o = 0, mean_squared_unc_jm_est_LT_o = 0, true_LT_o, con_jm_est_LT_o,mean_con_jm_est_LT_o = 0, mean_squared_con_jm_est_LT_o = 0;

    double unc_jm_est_HT_v, mean_unc_jm_est_HT_v = 0, mean_squared_unc_jm_est_HT_v = 0, true_HT_v, con_jm_est_HT_v, mean_con_jm_est_HT_v = 0, mean_squared_con_jm_est_HT_v = 0;
   double unc_jm_est_LT_v, mean_unc_jm_est_LT_v = 0, mean_squared_unc_jm_est_LT_v = 0, true_LT_v, con_jm_est_LT_v,mean_con_jm_est_LT_v = 0, mean_squared_con_jm_est_LT_v = 0;

   double se_temp, t_ratio_temp;

   int prop_rej_unc_jm_est_HT_o = 0,   prop_rej_con_jm_est_HT_o = 0;
   int prop_rej_unc_jm_est_LT_o = 0,   prop_rej_con_jm_est_L_o = 0;

   int prop_rej_unc_jm_est_HT_v = 0,   prop_rej_con_jm_est_HT_v = 0;
   int prop_rej_unc_jm_est_LT_v = 0,   prop_rej_con_jm_est_L_v = 0;

   double pval_rej_unc_jm_est_HT_o = 0,   pval_rej_con_jm_est_HT_o = 0;
   double pval_rej_unc_jm_est_LT_o = 0,   pval_rej_con_jm_est_L_o = 0;

   double pval_rej_unc_jm_est_HT_v = 0,   pval_rej_con_jm_est_HT_v = 0;
   double pval_rej_unc_jm_est_LT_v = 0,   pval_rej_con_jm_est_L_v = 0;


   true_HT_o = fix_ef_glmm_true(2);
   true_LT_o = fix_ef_glmm_true(3);

   true_HT_v = fix_ef_lme_true(2);
   true_LT_v = fix_ef_lme_true(3);

   double q;
   int simu_done = 0, simu_skipped = 0;
   do{
    Out<<"\n\n\n*************************************** simulation No."<<(simu_done+1)<<"\n";
   std::cout<<"simu_done="<<simu_done<<"\n";
    jm_gen_simu_data(rd,  nsub,  ni_vec, para_vec_true, cov_ranef_L_true,  dat);
    p_nsub = &nsub3;  p_ni_vec = &ni_vec;  pDat_ = &dat;
  // out_dat<<dat;


  //*************** separate modelling (GLMM and LME)
   // first compute unconstrained MLE, if already satisfy constraint, no need to compute constrained MLE.
   /*
   para_vec_glmm_ord = para_vec_glmm_ord_true; ConvFail = 0;
   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_glmm_d1_ord,  para_vec_glmm_ord, 11 , para_vec_glmm_lb,  para_vec_glmm_ub,  1,  1e-5,  20000, ConvFail );
   if(ConvFail==1) {  std::cout<<"skip since GLMM non-converge\n";  simu_skipped++; std::cout<<"simu_skipped = "<<simu_skipped<<"\n"; continue;}
  // Out<<"\n===GLMM, ConvFail = "<<ConvFail<<"\n";
 //  Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"unconstrained MLE of beta_o, HT, LT \n"<<trans(subm(para_vec_glmm_ord,range(2,3),range(0,0)));
 //  Out<<"trans(gradcdif(neg_loglike_glmm_d1_ord,para_vec_glmm_ord))\n"<<trans(gradcdif(neg_loglike_glmm_d1_ord,para_vec_glmm_ord));
   glmm_hess = HESScdif(neg_loglike_glmm_d1_ord,para_vec_glmm_ord);
   cov_glmm_mle_unc = inv(glmm_hess);
   cov_beta_o = subm(cov_glmm_mle_unc,range(2,3),range(2,3));
   t_ratio_unc_beta_o_HTi_min_LTi = para_vec_glmm_ord(2)/std::sqrt(cov_glmm_mle_unc(2,2));
   t_ratio_unc_beta_o_LTi = para_vec_glmm_ord(3)/std::sqrt(cov_glmm_mle_unc(3,3));
 // Out<<"cov_glmm_mle_unc\n"<<cov_glmm_mle_unc<<"cov_beta_o\n"<<cov_beta_o;
    for(int s1=2;s1<4;++s1)
  Out<<"sd_para_vec("<<s1<<") = "<<std::sqrt(cov_glmm_mle_unc(s1,s1))<<"\n";
   unc_mle_beta_o(0) = para_vec_glmm_ord(2);
   unc_mle_beta_o(1) = para_vec_glmm_ord(3);
   glmm_T_wald_ts = trans(unc_mle_beta_o)* inv(cov_beta_o) *unc_mle_beta_o ;
   pva_glmm_T_wald_ts = 1-scythe::pchisq(glmm_T_wald_ts,2);
 //  Out<<"glmm_T_wald_ts = "<<glmm_T_wald_ts<<"; pva_glmm_T_wald_ts = "<<pva_glmm_T_wald_ts<<"\n";
   if(pva_glmm_T_wald_ts<0.05) prop_rec_glmm_T_wald_ts++;


   // Out<<"unc_mle_beta_o\n"<<trans(unc_mle_beta_o);
  // Out<<"unc_mle_beta_o(0)<1e-5 & unc_mle_beta_o(1)<1e-5 = "<<(unc_mle_beta_o(0)<1e-5 & unc_mle_beta_o(1)<1e-5)<<"\n";

   if(unc_mle_beta_o(0)<1e-5 & unc_mle_beta_o(1)<1e-5) // unconstrained MLE is already within constraints
   {
    glmm_T_wald_os = glmm_T_wald_ts; // even thought test statistic value equal, p-val will differ since two-sided test is Chi-sq, One-sided test is Chi-bar-square
   }
   else
   {
      para_vec_glmm_ord = para_vec_glmm_ord_true; ConvFail = 0;
   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_glmm_d1_ord,  para_vec_glmm_ord, 11 , para_vec_glmm_ord_lb,  para_vec_glmm_ord_ub,  1,  1e-5,  20000, ConvFail );
   // Out<<"Constrained MLE trans(para_vec_glmm_ord)\n"<<trans(para_vec_glmm_ord);

   if(ConvFail==1) {  std::cout<<"skip since GLMM non-converge for inequality constrained mle\n";  simu_skipped++; std::cout<<"simu_skipped = "<<simu_skipped<<"\n"; continue;}
   con_mle_beta_o(0) = para_vec_glmm_ord(2);
   con_mle_beta_o(1) = para_vec_glmm_ord(3);
   glmm_T_wald_os = trans(con_mle_beta_o)* inv(cov_beta_o) *con_mle_beta_o ;
   }
   // one-sided wald test pval,
    q = std::acos(cov_beta_o(0,1)/sqrt(cov_beta_o(0,0)*cov_beta_o(1,1)) )/(2*pi);
    pva_glmm_T_wald_os = 1-q-0.5*scythe::pchisq(glmm_T_wald_os,1)-(0.5-q)*scythe::pchisq(glmm_T_wald_os,2);


   // LME model, separate analysis
   para_vec_lme_ord = para_vec_lme_ord_true; ConvFail = 0;
   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_lme_ana_ord,  para_vec_lme_ord, 2*6+1 , para_vec_lme_lb,  para_vec_lme_ub,  1,  1e-5,  20000, ConvFail );

   if(ConvFail==1) {  std::cout<<"skip since LME non-converge\n";  simu_skipped++; std::cout<<"simu_skipped = "<<simu_skipped<<"\n"; continue;}
 //  Out<<"LME, ConvFail = "<<ConvFail<<"\n";
 //  Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
    Out<<"unconstrained MLE beta_v, HT, LT\n"<<trans(subm(para_vec_lme_ord,range(2,3),range(0,0)));
      // Out<<"trans(gradcdif(neg_loglike_lme_d1,para_vec_lme_ord))\n"<<trans(gradcdif(neg_loglike_lme_ana_ord,para_vec_lme_ord));
   lme_hess = HESScdif(neg_loglike_lme_ana_ord,para_vec_lme_ord);
   cov_lme_mle_unc = inv(lme_hess);
   cov_beta_v = subm(cov_lme_mle_unc,range(2,3),range(2,3));
   t_ratio_unc_beta_v_HTi_min_LTi = para_vec_lme_ord(2)/std::sqrt(cov_lme_mle_unc(2,2));
   t_ratio_unc_beta_v_LTi = para_vec_lme_ord(3)/std::sqrt(cov_lme_mle_unc(3,3));
  // Out<<"cov_lme_mle_unc\n"<<cov_lme_mle_unc<<"cov_beta_v\n"<<cov_beta_v;
 for(int s1=2;s1<4;++s1)
  Out<<"sd_para_vec("<<s1<<") = "<<std::sqrt(cov_lme_mle_unc(s1,s1))<<"\n";
   unc_mle_beta_v(0) = para_vec_lme_ord(2);
   unc_mle_beta_v(1) = para_vec_lme_ord(3);
   lme_T_wald_ts = trans(unc_mle_beta_v)* inv(cov_beta_v) *unc_mle_beta_v ;
   pva_lme_T_wald_ts = 1-scythe::pchisq(lme_T_wald_ts,2);
 //  Out<<"lme_T_wald_ts = "<<lme_T_wald_ts<<"; pva_lme_T_wald_ts = "<<pva_lme_T_wald_ts<<"\n";
   if(pva_lme_T_wald_ts<0.05) prop_rec_lme_T_wald_ts++;

  //Out<<"unc_mle_beta_v\n"<<trans(unc_mle_beta_v);
 // Out<<"unc_mle_beta_v(0)<1e-5 & unc_mle_beta_v(1)<1e-5 = "<<(unc_mle_beta_v(0)<1e-5 & unc_mle_beta_v(1)<1e-5)<<"\n";

   if(unc_mle_beta_v(0)<1e-5 & unc_mle_beta_v(1)<1e-5) // unconstrained MLE is already within constraints
   {
    lme_T_wald_os = lme_T_wald_ts; // even thought test statistic value equal, p-val will differ since two-sided test is Chi-sq, One-sided test is Chi-bar-square
   }
   else
   {
   para_vec_lme_ord = para_vec_lme_ord_true; ConvFail = 0;
   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_lme_ana_ord,  para_vec_lme_ord, 2*6+1 , para_vec_lme_ord_lb,  para_vec_lme_ord_ub,  1,  1e-5,  20000, ConvFail );
   if(ConvFail==1) {  std::cout<<"skip since LME non-converge for inequality constrained mle\n";  simu_skipped++; std::cout<<"simu_skipped = "<<simu_skipped<<"\n"; continue;}
   //  Out<<"Constrained MLE trans(para_vec_lme_ord)\n"<<trans(para_vec_lme_ord);

   con_mle_beta_v(0) = para_vec_lme_ord(2);
   con_mle_beta_v(1) = para_vec_lme_ord(3);
   lme_T_wald_os = trans(con_mle_beta_v)* inv(cov_beta_v) *con_mle_beta_v ;
   }
      // one-sided wald test pval,
    q = std::acos(cov_beta_v(0,1)/sqrt(cov_beta_v(0,0)*cov_beta_v(1,1)) )/(2*pi);
    pva_lme_T_wald_os = 1-q-0.5*scythe::pchisq(lme_T_wald_os,1)-(0.5-q)*scythe::pchisq(lme_T_wald_os,2);
    // Out<<"lme_T_wald_os = "<<lme_T_wald_os<<"; pva_lme_T_wald_os = "<<pva_lme_T_wald_os<<"\n";

    out_log<<con_mle_beta_o(0)<<","<<unc_mle_beta_o(0)<<","<<t_ratio_unc_beta_o_HTi_min_LTi<<","<<con_mle_beta_o(1)<<","<<unc_mle_beta_o(1)<<","<<t_ratio_unc_beta_o_LTi<<","<<glmm_T_wald_ts<<","<<glmm_T_wald_os<<","<<pva_glmm_T_wald_ts<<","<<pva_glmm_T_wald_os<<","<<con_mle_beta_v(0)<<","<<unc_mle_beta_v(0)<<","<<t_ratio_unc_beta_v_HTi_min_LTi<<","<<con_mle_beta_v(1)<<","<<unc_mle_beta_v(1)<<","<<t_ratio_unc_beta_v_LTi<<","<<lme_T_wald_ts<<","<<lme_T_wald_os<<","<<pva_lme_T_wald_ts<<","<<pva_lme_T_wald_os<<"\n";
   */

    // joint model
    para_vec_simu_ord = para_vec_true_ord;
    //Out<<"para_vec_simu_ord\n"<<trans(para_vec_simu_ord);
   p_GH_points = &aQHPts9;
    //Out<<"initial neg_loglike_lme_ana_glmm_GH_ord(para_vec_simu_ord)\n"<<neg_loglike_lme_ana_glmm_GH_ord(para_vec_simu_ord)<<"\n";
     ConvFail = 0;
   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_lme_ana_glmm_GH_ord,  para_vec_simu_ord, Two_num_para_plus1 , lb,  ub,  1,  bobyqa_stop_criteria,  20000, ConvFail );

   if(ConvFail==1) {  std::cout<<"skip since joint model non-converge\n";  simu_skipped++; std::cout<<"simu_skipped = "<<simu_skipped<<"\n"; continue;}
 //  Out<<"LME, ConvFail = "<<ConvFail<<"\n";
 //  Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
 unc_mle_HT_LT(0) = para_vec_simu_ord(2);
 unc_mle_HT_LT(1) = para_vec_simu_ord(3);
 unc_mle_HT_LT(2) = para_vec_simu_ord(6);
 unc_mle_HT_LT(3) = para_vec_simu_ord(7);

 unc_jm_est_HT_o = para_vec_simu_ord(2) + para_vec_simu_ord(3);
 mean_unc_jm_est_HT_o = mean_unc_jm_est_HT_o + unc_jm_est_HT_o;
 mean_squared_unc_jm_est_HT_o = mean_squared_unc_jm_est_HT_o + unc_jm_est_HT_o*unc_jm_est_HT_o;
   unc_jm_est_LT_o =   para_vec_simu_ord(3);
 mean_unc_jm_est_LT_o = mean_unc_jm_est_LT_o + unc_jm_est_LT_o;
 mean_squared_unc_jm_est_LT_o = mean_squared_unc_jm_est_LT_o + unc_jm_est_LT_o*unc_jm_est_LT_o;

 unc_jm_est_HT_v = para_vec_simu_ord(6) + para_vec_simu_ord(7);
 mean_unc_jm_est_HT_v = mean_unc_jm_est_HT_v + unc_jm_est_HT_v;
 mean_squared_unc_jm_est_HT_v = mean_squared_unc_jm_est_HT_v + unc_jm_est_HT_v*unc_jm_est_HT_v;
   unc_jm_est_LT_v =   para_vec_simu_ord(7);
 mean_unc_jm_est_LT_v = mean_unc_jm_est_LT_v + unc_jm_est_LT_v;
 mean_squared_unc_jm_est_LT_v = mean_squared_unc_jm_est_LT_v + unc_jm_est_LT_v*unc_jm_est_LT_v;

   Out<<"JM, unc_mle_HT_LT\n"<<trans(unc_mle_HT_LT);
   Out<<"trans(gradcdif(neg_loglike_lme_d1,para_vec_lme_ord))\n"<<trans(gradcdif(neg_loglike_lme_ana_glmm_GH_ord,para_vec_simu_ord));
   jm_hess = HESScdif(neg_loglike_lme_ana_glmm_GH_ord,para_vec_simu_ord);
   cov_jm_mle_unc = inv(jm_hess);
   set_subm(cov_mle_HT_LT,range(0,1),range(0,1)) = subm(cov_jm_mle_unc,range(2,3),range(2,3));
   set_subm(cov_mle_HT_LT,range(2,3),range(2,3)) = subm(cov_jm_mle_unc,range(6,7),range(6,7));

   // Out<<"cov_mle_HT_LT\n"<<cov_mle_HT_LT;
   for(int s1=0;s1<4;++s1)
  Out<<"sd_mle_HT_LT("<<s1<<") = "<<std::sqrt(cov_mle_HT_LT(s1,s1))<<"\n";
  inv_cov_mle_HT_LT = inv(cov_mle_HT_LT);
   jm_T_wald_ts = trans(unc_mle_HT_LT)* inv_cov_mle_HT_LT*unc_mle_HT_LT ;
   pva_jm_T_wald_ts = 1-scythe::pchisq(jm_T_wald_ts,4);
   if(pva_jm_T_wald_ts<0.05) prop_rec_jm_T_wald_ts++;
   Out<<"pva_jm_T_wald_ts = "<<pva_jm_T_wald_ts<<"; jm_T_wald_ts = "<<jm_T_wald_ts<<"\n";





    //Out<<"para_vec_true_ord\n"<<trans(para_vec_true_ord);
    //Out<<"lb\n"<<trans(lb);
    //Out<<"ub_ord\n"<<trans(ub_ord);

  // std::cout<<" unc_mle_HT_LT(0)<1e-5 & unc_mle_HT_LT(1)<1e-5 & unc_mle_HT_LT(2)<1e-5 & unc_mle_HT_LT(3)<1e-5 = "<<(unc_mle_HT_LT(0)<1e-5 & unc_mle_HT_LT(1)<1e-5 & unc_mle_HT_LT(2)<1e-5 & unc_mle_HT_LT(3)<1e-5)<<"\n";


    if(unc_mle_HT_LT(0)<1e-5 & unc_mle_HT_LT(1)<1e-5 & unc_mle_HT_LT(2)<1e-5 & unc_mle_HT_LT(3)<1e-5 ) // unconstrained MLE is already within constraints
   {
    jm_T_wald_os = jm_T_wald_ts; // even thought test statistic value equal, p-val will differ since two-sided test is Chi-sq, One-sided test is Chi-bar-square
    con_jm_est_HT_o = para_vec_simu_ord(2) + para_vec_simu_ord(3);
 mean_con_jm_est_HT_o = mean_con_jm_est_HT_o + con_jm_est_HT_o;
 mean_squared_con_jm_est_HT_o = mean_squared_con_jm_est_HT_o + con_jm_est_HT_o*con_jm_est_HT_o;
   con_jm_est_LT_o = para_vec_simu_ord(3);
 mean_con_jm_est_LT_o = mean_con_jm_est_LT_o + con_jm_est_LT_o;
 mean_squared_unc_jm_est_LT_o = mean_squared_unc_jm_est_LT_o + con_jm_est_LT_o*con_jm_est_LT_o;

   con_jm_est_HT_v = para_vec_simu_ord(6) + para_vec_simu_ord(7);
 mean_con_jm_est_HT_v = mean_con_jm_est_HT_v + con_jm_est_HT_v;
 mean_squared_con_jm_est_HT_v = mean_squared_con_jm_est_HT_v + con_jm_est_HT_v*unc_jm_est_HT_v;
   con_jm_est_LT_v =   para_vec_simu_ord(7);
 mean_con_jm_est_LT_v = mean_con_jm_est_LT_v + con_jm_est_LT_v;
 mean_squared_con_jm_est_LT_v = mean_squared_con_jm_est_LT_v + con_jm_est_LT_v*unc_jm_est_LT_v;

   }
   else
   {
      para_vec_simu_ord = para_vec_true_ord; ConvFail = 0;
  //std::cout<<"before constrained MLE, para_vec_simu_ord\n"<<trans(para_vec_simu_ord);

   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_lme_ana_glmm_GH_ord,  para_vec_simu_ord, Two_num_para_plus1 , lb,  ub_ord,  1,  bobyqa_stop_criteria,  20000, ConvFail );

   if(ConvFail==1) {  std::cout<<"skip since JM has non-converge for inequality constrained mle\n";  simu_skipped++; std::cout<<"simu_skipped = "<<simu_skipped<<"\n"; continue;}

    con_mle_HT_LT(0) = para_vec_simu_ord(2);
 con_mle_HT_LT(1) = para_vec_simu_ord(3);
 con_mle_HT_LT(2) = para_vec_simu_ord(6);
 con_mle_HT_LT(3) = para_vec_simu_ord(7);
   Out<<"JM, con_mle_HT_LT\n"<<trans(con_mle_HT_LT);


      con_jm_est_HT_o = para_vec_simu_ord(2) + para_vec_simu_ord(3);
 mean_con_jm_est_HT_o = mean_con_jm_est_HT_o + con_jm_est_HT_o;
 mean_squared_con_jm_est_HT_o = mean_squared_con_jm_est_HT_o + con_jm_est_HT_o*con_jm_est_HT_o;
   con_jm_est_LT_o = para_vec_simu_ord(3);
 mean_con_jm_est_LT_o = mean_con_jm_est_LT_o + con_jm_est_LT_o;
 mean_squared_unc_jm_est_LT_o = mean_squared_unc_jm_est_LT_o + con_jm_est_LT_o*con_jm_est_LT_o;

   con_jm_est_HT_v = para_vec_simu_ord(6) + para_vec_simu_ord(7);
 mean_con_jm_est_HT_v = mean_con_jm_est_HT_v + con_jm_est_HT_v;
 mean_squared_con_jm_est_HT_v = mean_squared_con_jm_est_HT_v + con_jm_est_HT_v*unc_jm_est_HT_v;
   con_jm_est_LT_v =   para_vec_simu_ord(7);
 mean_con_jm_est_LT_v = mean_con_jm_est_LT_v + con_jm_est_LT_v;
 mean_squared_con_jm_est_LT_v = mean_squared_con_jm_est_LT_v + con_jm_est_LT_v*unc_jm_est_LT_v;


   jm_T_wald_os = trans(con_mle_HT_LT)* inv_cov_mle_HT_LT *con_mle_HT_LT;
   }
   // one-sided wald test pval,
   L_cov_mle_chi_bar_cdf_simu = chol(cov_mle_HT_LT);
   chi_bar_d4_cdf_weights(cov_mle_HT_LT, inv_cov_mle_HT_LT,rd,num_chi_bar_sam, proj_of_mle_chi_bar_cdf_simu,  L_cov_mle_chi_bar_cdf_simu, std_mvn_sam, L_mvn_sam, lb_chi_bar_cdf_simu, ub_chi_bar_cdf_simu, chi_bar_weights);
   //  std::cout<<"chi_bar_weights\n"<<trans(chi_bar_weights);
   pva_jm_T_wald_os = 1-chi_bar_d4_cdf(jm_T_wald_os,chi_bar_weights);

   if(pva_jm_T_wald_os<0.05) prop_rec_jm_T_wald_os++;
   Out<<"pva_jm_T_wald_os = "<<pva_jm_T_wald_os<<"; jm_T_wald_os = "<<jm_T_wald_os<<"\n";


   for(s1_para=0;s1_para<4;++s1_para)
   {
   se_unc_mle_HT_LT(s1_para) = std::sqrt(cov_mle_HT_LT(s1_para,s1_para));
   t_ratio_unc_mle_HT_LT(s1_para) = unc_mle_HT_LT(s1_para)/se_unc_mle_HT_LT(s1_para);
   t_ratio_con_mle_HT_LT(s1_para) = con_mle_HT_LT(s1_para)/se_unc_mle_HT_LT(s1_para);
   uni_p_val_unc_mle_HT_LT(s1_para) = 2*(1-scythe::pnorm(std::abs( t_ratio_unc_mle_HT_LT(s1_para)),0,1));
   uni_p_val_con_mle_HT_LT(s1_para) = scythe::pnorm(t_ratio_con_mle_HT_LT(s1_para),0,1);
   }

     out_log<<unc_jm_est_HT_v<<","<<unc_jm_est_LT_v<<","<<unc_jm_est_HT_o<<","<<unc_jm_est_LT_o<<","<<con_jm_est_HT_v<<","<<con_jm_est_LT_v<<","<<con_jm_est_HT_o<<","<<con_jm_est_LT_o<<","<<pva_jm_T_wald_os<<","<<pva_jm_T_wald_ts<<","<<uni_p_val_unc_mle_HT_LT(0)<<","<<uni_p_val_unc_mle_HT_LT(1)<<","<<uni_p_val_unc_mle_HT_LT(2)<<","<<uni_p_val_unc_mle_HT_LT(3)<<","<<uni_p_val_con_mle_HT_LT(0)<<","<<uni_p_val_con_mle_HT_LT(1)<<","<<uni_p_val_con_mle_HT_LT(2)<<","<<uni_p_val_con_mle_HT_LT(3)<<"\n";

  simu_done++;
   } while(simu_done < num_simu_replication);

   Out<<"\n\n\nPercentage of H_0 rej, Two-sided Wald, prop_rec_jm_T_wald_ts\n";
   Out<<( 100*prop_rec_jm_T_wald_ts/double(simu_done) )<<"\n\n";
   Out<<"Percentage of H_0 rej, One-sided Wald, prop_rec_jm_T_wald_os\n";
   Out<<( 100*prop_rec_jm_T_wald_os/double(simu_done) )<<"\n\n";

   /*
   mean_unc_jm_est_HT_o = mean_unc_jm_est_HT_o/double(simu_done);
   mean_unc_jm_est_HT_v = mean_unc_jm_est_HT_v/double(simu_done);
  mean_squared_unc_jm_est_HT_o = mean_squared_unc_jm_est_HT_o/double(simu_done);
  mean_squared_unc_jm_est_HT_v = mean_squared_unc_jm_est_HT_v/double(simu_done);

     mean_con_jm_est_HT_o = mean_con_jm_est_HT_o/double(simu_done);
   mean_con_jm_est_HT_v = mean_con_jm_est_HT_v/double(simu_done);
  mean_squared_con_jm_est_HT_o = mean_squared_con_jm_est_HT_o/double(simu_done);
  mean_squared_con_jm_est_HT_v = mean_squared_con_jm_est_HT_v/double(simu_done);

   Out<<"\n**Bias Unc HT_o = "<<(mean_unc_jm_est_HT_o - true_HT_o)<<"\n\n";
   Out<<"SE Unc HT_o = "<<std::sqrt(mean_squared_unc_jm_est_HT_o - mean_unc_jm_est_HT_o*mean_unc_jm_est_HT_o)<<"\n\n";
   Out<<"MSE Unc HT_o = "<<( (mean_unc_jm_est_HT_o - true_HT_o)*(mean_unc_jm_est_HT_o - true_HT_o) +(mean_squared_unc_jm_est_HT_o - mean_unc_jm_est_HT_o*mean_unc_jm_est_HT_o) )<<"\n\n";

   Out<<"\n**Bias Con HT_o = "<<(mean_con_jm_est_HT_o - true_HT_o)<<"\n\n";
   Out<<"SE Con HT_o = "<<std::sqrt(mean_squared_con_jm_est_HT_o - mean_con_jm_est_HT_o*mean_con_jm_est_HT_o)<<"\n\n";
   Out<<"MSE Con HT_o = "<<( (mean_con_jm_est_HT_o - true_HT_o)*(mean_con_jm_est_HT_o - true_HT_o)+(mean_squared_con_jm_est_HT_o - mean_con_jm_est_HT_o*mean_con_jm_est_HT_o) )<<"\n\n";

      Out<<"\n**Bias Unc LT_o = "<<(mean_unc_jm_est_LT_o - true_LT_o)<<"\n\n";
   Out<<"SE Unc LT_o = "<<std::sqrt(mean_squared_unc_jm_est_LT_o - mean_unc_jm_est_LT_o*mean_unc_jm_est_LT_o)<<"\n\n";
   Out<<"MSE Unc LT_o = "<<( (mean_unc_jm_est_LT_o - true_LT_o)*(mean_unc_jm_est_LT_o - true_LT_o)+(mean_squared_unc_jm_est_LT_o - mean_unc_jm_est_LT_o*mean_unc_jm_est_LT_o) )<<"\n\n";

   Out<<"\n**Bias Con LT_o = "<<(mean_con_jm_est_LT_o - true_LT_o)<<"\n\n";
   Out<<"SE Con LT_o = "<<std::sqrt(mean_squared_con_jm_est_LT_o - mean_con_jm_est_LT_o*mean_con_jm_est_LT_o)<<"\n\n";
   Out<<"MSE Con LT_o = "<<( (mean_con_jm_est_LT_o - true_LT_o)*(mean_con_jm_est_LT_o - true_LT_o)+(mean_squared_con_jm_est_LT_o - mean_con_jm_est_LT_o*mean_con_jm_est_LT_o) )<<"\n\n";

   Out<<"\n**Bias Unc HT_v = "<<(mean_unc_jm_est_HT_v - true_HT_v)<<"\n\n";
   Out<<"SE Unc HT_v = "<<std::sqrt(mean_squared_unc_jm_est_HT_v - mean_unc_jm_est_HT_v*mean_unc_jm_est_HT_v)<<"\n\n";
   Out<<"MSE Unc HT_v = "<<( (mean_unc_jm_est_HT_v - true_HT_v)*(mean_unc_jm_est_HT_v - true_HT_v) +(mean_squared_unc_jm_est_HT_v - mean_unc_jm_est_HT_v*mean_unc_jm_est_HT_v) )<<"\n\n";

   Out<<"\n**Bias Con HT_v = "<<(mean_con_jm_est_HT_v - true_HT_v)<<"\n\n";
   Out<<"SE Con HT_v = "<<std::sqrt(mean_squared_con_jm_est_HT_v - mean_con_jm_est_HT_v*mean_con_jm_est_HT_v)<<"\n\n";
   Out<<"MSE Con HT_v = "<<( (mean_con_jm_est_HT_v - true_HT_v)*(mean_con_jm_est_HT_v - true_HT_v)+(mean_squared_con_jm_est_HT_v - mean_con_jm_est_HT_v*mean_con_jm_est_HT_v) )<<"\n\n";

      Out<<"\n**Bias Unc LT_v = "<<(mean_unc_jm_est_LT_v - true_LT_v)<<"\n\n";
   Out<<"SE Unc LT_v = "<<std::sqrt(mean_squared_unc_jm_est_LT_v - mean_unc_jm_est_LT_v*mean_unc_jm_est_LT_v)<<"\n\n";
   Out<<"MSE Unc LT_v = "<<( (mean_unc_jm_est_LT_v - true_LT_v)*(mean_unc_jm_est_LT_v - true_LT_v)+(mean_squared_unc_jm_est_LT_v - mean_unc_jm_est_LT_v*mean_unc_jm_est_LT_v) )<<"\n\n";

   Out<<"\n**Bias Con LT_v = "<<(mean_con_jm_est_LT_v - true_LT_v)<<"\n\n";
   Out<<"SE Con LT_v = "<<std::sqrt(mean_squared_con_jm_est_LT_v - mean_con_jm_est_LT_v*mean_con_jm_est_LT_v)<<"\n\n";
   Out<<"MSE Con LT_v = "<<( (mean_con_jm_est_LT_v - true_LT_v)*(mean_con_jm_est_LT_v - true_LT_v)+(mean_squared_con_jm_est_LT_v - mean_con_jm_est_LT_v*mean_con_jm_est_LT_v) )<<"\n\n";
 */

  return 0;
}

  //  p_GH_points = &aQHPts13; para_vec_lme = para_vec_lme_true; ConvFail = 0;
 //  Out<<"\n\n GH13, initial value GH_neg_loglike_lme_d1(para_vec_lme) = "<< GH_neg_loglike_lme_d1(para_vec_lme)<<"\n";
 //  Out<<"para_vec_lme\n"<<trans(para_vec_lme);
 //     p_GH_points = &aQHPts19;
 //  Out<<"GH19, initial value GH_neg_loglike_lme_d1(para_vec_lme) = "<< GH_neg_loglike_lme_d1(para_vec_lme)<<"\n";
  //  Out<<"neg_loglike_lme_ana(para_vec_lme) = "<<neg_loglike_lme_ana(para_vec_lme)<<"\n";





    /*
   min_neg_loglike_unc = find_min_bobyqaConvgFail(neg_loglike_lme_d1,  para_vec_lme, 2*6+1 , para_vec_lme_lb,  para_vec_lme_ub,  1,  1e-4,  200, ConvFail );


   Out<<"LME, ConvFail = "<<ConvFail<<"\n";
   Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"trans(para_vec_lme)\n"<<trans(para_vec_lme);
   Out<<"trans(gradcdif(neg_loglike_lme_d1,para_vec_lme))\n"<<trans(gradcdif(neg_loglike_lme_d1,para_vec_lme));
   lme_hess = HESScdif(neg_loglike_lme_d1,para_vec_lme);
   cov_lme_mle_unc = inv(lme_hess);
   for(int s1=0;s1<4;++s1)
   Out<<"sd_para_vec("<<s1<<") = "<<std::sqrt(cov_lme_mle_unc(s1,s1))<<"\n";



   para_vec_lme_Lap = para_vec_lme_true_Lap; ConvFail = 0;
   Out<<"para_vec_lme_Lap\n"<<trans(para_vec_lme_Lap);
   Out<<"\n\n initial value Lap_neg_loglike_lme_d1(para_vec_lme_Lap) = "<< Lap_neg_loglike_lme_d1(para_vec_lme_Lap)<<"\n";
   min_neg_loglike_unc = find_min_bobyqaConvgFail(Lap_neg_loglike_lme_d1,  para_vec_lme_Lap, 2*6+1 , para_vec_lme_lb,  para_vec_lme_ub, 1,  1e-4,  200, ConvFail );
   Out<<"LME Laplace, ConvFail = "<<ConvFail<<"\n";
   Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"trans(para_vec_lme)\n"<<trans(para_vec_lme);
   Out<<"trans(gradcdif(neg_loglike_lme_d1,para_vec_lme))\n"<<trans(gradcdif(Lap_neg_loglike_lme_d1,para_vec_lme));
   lme_hess = HESScdif(Lap_neg_loglike_lme_d1,para_vec_lme);
   cov_lme_mle_unc = inv(lme_hess);
   for(int s1=0;s1<4;++s1)
   Out<<"sd_para_vec("<<s1<<") = "<<std::sqrt(cov_lme_mle_unc(s1,s1))<<"\n";


 //infLoglike_para_vec = 0    ,     0 ,        0 ,        0 ,       -1   ,      0  ,       0   ,      0, -3.50655789732, -0.673536823983 ,        0, -3.10730404921;
 //  para_vec_simu = infLoglike_para_vec;


    p_rd = &rd;  dbcolvec std_mvn_sam(dim_ranf), ranf_mvn_sam(dim_ranf); p_std_mvn_sam = &std_mvn_sam; p_ranf_mvn_sam = &ranf_mvn_sam;
 int num_mc = 3000; p_num_mc = &num_mc;
 dbcolvec ranef_sam_mean(dim_ranf); p_sam_mean = &ranef_sam_mean;
  dbmat ranef_sam_cov(dim_ranf,dim_ranf); p_sam_cov = &ranef_sam_cov;
  double neg_loglike_U95; p_neg_loglike_U95 = &neg_loglike_U95;
  double neg_loglike_L95; p_neg_loglike_L95 = &neg_loglike_L95;
  Out<<"\n\n*** MC_neg_loglike(para_vec_true) = "<<MC_neg_loglike(para_vec_true)<<"\n";
 Out<<"MC L95 = "<<neg_loglike_L95<<"\n";
 Out<<"MC U95 = "<<neg_loglike_U95<<"\n";
  Out<<"*p_sam_mean = "<<*p_sam_mean;
  Out<<"*p_sam_cov = "<<*p_sam_cov<<"\n";


 // p_GH_points = &aQHPts5;  Out<<"GH5: "<<GH_neg_loglike_direct_d2(para_vec_true)<<"\n";
  //  p_GH_points = &aQHPts7; Out<<"GH7: "<<GH_neg_loglike_direct_d2(para_vec_true)<<"\n";
  // p_GH_points = &aQHPts9; Out<<"GH9: "<<GH_neg_loglike_direct_d2(para_vec_true)<<"\n";
  //
  p_GH_points = &aQHPts19; Out<<"GH19: "<<GH_neg_loglike_direct_d2(para_vec_true)<<"\n";



   para_vec_simu_Lap = para_vec_true_Lap; ConvFail = 0;
  Out<<"trans(para_vec_simu_Lap initial)\n"<<trans(para_vec_simu_Lap);
 Out<<"Lap_neg_loglike(para_vec_Lap initial): "<<Lap_neg_loglike(para_vec_simu_Lap)<<"\n";
 min_neg_loglike_unc =  find_min_bobyqaConvgFail(Lap_neg_loglike,  para_vec_simu_Lap, Two_num_para_plus1 , lb,  ub,  1.0,  1e-4,  200, ConvFail );
   Out<<"ConvFail = "<<ConvFail<<"\n";
   grad_at_unc_mle = gradcdif(Lap_neg_loglike,para_vec_simu_Lap);
   Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"trans(para_vec_simu_Lap)\n"<<trans(para_vec_simu_Lap);
   Out<<"trans(grad_at_unc_mle)\n"<<trans(grad_at_unc_mle);

    para_vec_simu_Lap = para_vec_true_Lap; ConvFail = 0;
  min_neg_loglike_unc =  find_minConvFail(bfgs_search_strategy(), objective_delta_stop_strategy(1e-7), Lap_neg_loglike, Grad_Lap_neg_loglike, para_vec_simu_Lap, -BIGNUM,ConvFail,100);
   Out<<"\n\n\nConvFail = "<<ConvFail<<"\n";
   grad_at_unc_mle = gradcdif(Lap_neg_loglike,para_vec_simu_Lap);
   Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"trans(para_vec_simu_Lap)\n"<<trans(para_vec_simu_Lap);
   Out<<"trans(grad_at_unc_mle)\n"<<trans(grad_at_unc_mle);



  // para_vec_simu = 3;  set_subm(para_vec_simu,range(9,11),range(0,0)) = cov_ranef_Lvec_true;


   para_vec_simu = para_vec_true;
   p_GH_points = &aQHPts13; ConvFail = 0;
  Out<<"trans(para_vec_simu initial)\n"<<trans(para_vec_simu);
  Out<<"GH_neg_loglike_direct_d2(para_vec_simu) initial: "<<GH_neg_loglike_direct_d2(para_vec_simu)<<"\n";
 min_neg_loglike_unc =  find_min_bobyqaConvgFail(GH_neg_loglike_direct_d2,  para_vec_simu, Two_num_para_plus1 , lb,  ub,  0.3,  1e-4,  200, ConvFail );
   Out<<"ConvFail = "<<ConvFail<<"\n";
   grad_at_unc_mle = gradcdif(GH_neg_loglike_direct_d2,para_vec_simu);
   Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"trans(para_vec_simu)\n"<<trans(para_vec_simu);
   Out<<"trans(grad_at_unc_mle)\n"<<trans(grad_at_unc_mle);


    p_GH_points = &aQHPts13; ConvFail = 0; para_vec_simu = para_vec_true;
     min_neg_loglike_unc = -1.0*find_minConvFail(bfgs_search_strategy(), objective_delta_stop_strategy(1e-6).be_verbose(), GH_neg_loglike_direct_d2, Grad_GH_neg_loglike_direct_d2, para_vec_simu, -1e99,ConvFail,2000);  //
   Out<<"min_neg_loglike_unc: "<<min_neg_loglike_unc<<"\n";
   Out<<"trans(para_vec_simu)\n"<<trans(para_vec_simu);
   Out<<"trans(grad_at_unc_mle)\n"<<trans(grad_at_unc_mle);
    */



    /*
   p_GH_points = &aQHPts13; Out<<"GH neg loglike = "<<GH_neg_loglike_direct_d2(para_vec_true)<<"\n";
   p_rd = &rd;  dbcolvec std_mvn_sam(dim_ranf), ranf_mvn_sam(dim_ranf); p_std_mvn_sam = &std_mvn_sam; p_ranf_mvn_sam = &ranf_mvn_sam;
 int num_mc = 3000; p_num_mc = &num_mc;
 dbcolvec ranef_sam_mean(dim_ranf); p_sam_mean = &ranef_sam_mean;
  dbmat ranef_sam_cov(dim_ranf,dim_ranf); p_sam_cov = &ranef_sam_cov;
  double neg_loglike_U95; p_neg_loglike_U95 = &neg_loglike_U95;
  double neg_loglike_L95; p_neg_loglike_L95 = &neg_loglike_L95;
  Out<<"\n\n*** MC_neg_loglike(para_vec_true) = "<<MC_neg_loglike(para_vec_true)<<"\n";
 Out<<"MC L95 = "<<neg_loglike_L95<<"\n";
 Out<<"MC U95 = "<<neg_loglike_U95<<"\n";
  Out<<"*p_sam_mean = "<<*p_sam_mean;
  Out<<"*p_sam_cov = "<<*p_sam_cov<<"\n";


   Out<<"neg_loglike_glmm_d1_ord(para_vec_glmm_ord_true) = "<< neg_loglike_glmm_d1_ord(para_vec_glmm_ord_true)<<"\n";

    Out<<"neg_loglike_lme_ana_ord(para_vec_lme_ord) = "<<neg_loglike_lme_ana_ord(para_vec_lme_ord)<<"\n";
   para_vec_lme_ord = para_vec_lme_ord_true;

   */
