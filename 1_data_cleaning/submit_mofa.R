#calculate PEER factors using MOFA
load("./protein_flag2_all.RData")

#BiocManager::install("MOFA2")
library(MOFA2)
dat_input = t(protein_clean_final[,-1])
colnames(dat_input) = protein_clean_final$SampleId
dat_input = list(dat_input)
MOFAobject <- create_mofa(dat_input)
plot_data_overview(MOFAobject)
#define data options
data_opts <- get_default_data_options(MOFAobject)
#define model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors = 10
head(model_opts)
#define training options
train_opts <- get_default_training_options(MOFAobject)

#Build and train the MOFA object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
outfile = file.path("./protein","ARIC_proteomics.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject,outfile,use_basilisk = T)
