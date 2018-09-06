/* Set control parameters. */
LIKELIHOOD_FUNCTION_OUTPUT = 5;
MAXIMUM_ITERATIONS_PER_VARIABLE = 10000;

/* Initialize (pseudo)random-number generator. */
SetParameter(RANDOM_SEED, random_seed, 0);

/* Set up query compartment data subset. */
DataSet quer_data_subset = ReadDataFile(quer_seq_file);
DataSetFilter quer_data_subset_filt = CreateFilter(quer_data_subset, 1);
HarvestFrequencies(quer_pis, quer_data_subset, 1, 1, 1);

/* Set up reference compartment data subset. */
DataSet ref_data_subset = ReadDataFile(ref_seq_file);
DataSetFilter ref_data_subset_filt = CreateFilter(ref_data_subset, 1);
HarvestFrequencies(ref_pis, ref_data_subset, 1, 1, 1);

/* Set up model. */
global kappa_inverse;

/* Set up query compartment submodel. */
global f0;
f0 :< 1;
global f_aux;
f_aux :< 1;
site_class_freqs = {{f0, (1.0-f0)*f_aux, (1.0-f0)*(1.0-f_aux)*f0/(f0 + (1.0-f0)*f_aux), (1.0-f0)*(1.0-f_aux)*(1.0-f0)*f_aux/(f0 + (1.0-f0)*f_aux)}};
site_class_vals = {{0, 1, 2, 3}};
category site_class = (4, site_class_freqs, MEAN, , site_class_vals, 0, 3);
global zeta0;
zeta0 :< 1;
global zeta2;
zeta2 :> 1;
global zeta_bgrnd;
zeta_bgrnd := ((site_class == 0)+(site_class == 2))*zeta0 + ((site_class == 1)+(site_class == 3));
global zeta_fgrnd;
zeta_fgrnd := (site_class == 0)*zeta0 + (site_class == 1) + ((site_class == 2)+(site_class == 3))*zeta2;
quer_mat = {{*, kappa_inverse*t, t, kappa_inverse*t}
            {kappa_inverse*t, *, kappa_inverse*t, t}
            {t, kappa_inverse*t, *, kappa_inverse*t}
            {kappa_inverse*t, t, kappa_inverse*t, *}};
Model quer_submod = (quer_mat, quer_pis, 1);

/* Set up reference compartment submodel. */
ref_mat = {{*, kappa_inverse*t, t, kappa_inverse*t}
           {kappa_inverse*t, *, kappa_inverse*t, t}
           {t, kappa_inverse*t, *, kappa_inverse*t}
           {kappa_inverse*t, t, kappa_inverse*t, *}};
Model ref_submod = (ref_mat, ref_pis, 1);

/* Fit model to data sets. */
log_Ls = {fit_repl_count, 1};
for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1) {
  kappa_inverse = Random(0.0, 10.0);
  f0 = Random(0.0, 1.0);
  f_aux = Random(0.0, 1.0);
  zeta0 = Random(0.0, 1.0);
  zeta2 = Random(1.0, 10.0);
  UseModel(quer_submod);
  Tree quer_hyphy_tree = tree;
  UseModel(ref_submod);
  Tree ref_hyphy_tree = tree;
  ReplicateConstraint("this1.?.t := zeta_bgrnd*this2.?.t", quer_hyphy_tree, ref_hyphy_tree);
  ExecuteCommands("quer_hyphy_tree."+fgrnd_branch_name+".t := zeta_fgrnd*ref_hyphy_tree."+fgrnd_branch_name+".t;");
  LikelihoodFunction L = (quer_data_subset_filt, quer_hyphy_tree, ref_data_subset_filt, ref_hyphy_tree);
  Optimize(MLEs, L);
  log_L = MLEs[1][0];
  log_Ls[repl_i] = log_L;
  if (repl_i == 0 || best_log_L < log_L) {
    best_repl_i = repl_i;
    best_log_L = log_L;
    best_MLEs = MLEs;
  }
  fprintf(stdout, "REPL ", repl_i, "\n");
  fprintf(stdout, "log_L:", Format(log_L, 0, 16));
  fprintf(stdout, L, "\n\n");
}
sum = 0.0;
for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1) sum = sum+log_Ls[repl_i];
mean_log_L = sum/fit_repl_count;
if (fit_repl_count > 1) {
  sum = 0.0;
  for (repl_i = 0; repl_i < fit_repl_count; repl_i = repl_i+1) sum = sum + (log_Ls[repl_i]-mean_log_L)^2;
  stdev_log_L = Sqrt(sum/(fit_repl_count-1));
} else
  stdev_log_L = 0.0;

/* Write results. */
fprintf(res_file, "Results File: ", res_file, "\n");
fprintf(res_file, "BEST REPL: ", best_repl_i, "\n");
fprintf(res_file, "BEST LOG-L: ", Format(best_log_L, 0, 6), "\n");
