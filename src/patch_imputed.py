import os
import re
import pandas
import gzip
import pyarrow.parquet as pq
import subprocess
import shlex
import pickle
import numpy
import logging

logging.basicConfig(format='%(asctime)s %(message)s')
logging.getLogger().setLevel(logging.INFO)

def read_imputed_gwas(gwas_filename):
  _open = gzip.open if gwas_filename.endswith(".gz") else open
  with _open(gwas_filename) as gwas:
    df = pandas.read_csv(gwas, sep="\t")
  return df


##################################
##   FILLING WITH SAMPLE SIZE   ## 
##################################

# file with sample sizes (controls and cases when appropriate).
# possible file schema:
# PHENOTYPE_NAME_AS_IN_FILE\tQUANT\tCASES\tCONTROLS

def read_sample_size(sample_size_file, header=True):
  dd = {}
  with open(sample_size_file) as ssf:
    if header: 
      ssf.readline()
    for line in ssf:
      comps = line.strip().split()
      dd[comps[0]] = {'binary': bool(int(comps[1])), \
                      'sample_size': int(comps[3]), \
                      'cases': int(comps[2]) if bool(int(comps[1])) else None}
  return(dd)


def is_sample_size_available(gwas_df):
  # returns True if not all of the cells are NA's
  return(not gwas_df['sample_size'].isnull().all())


# calculate median when sample size is available
def get_median_sample_size(gwas_df):
  median_ss = gwas_df["sample_size"][~gwas_df["sample_size"].isnull()].median()
  return(median_ss)
 

def complete_missing_sample_sizes(gwas_df, phenotype, sample_size_d):
  if is_sample_size_available(gwas_df):
    median_sample_size = get_median_sample_size(gwas_df)
    print(median_sample_size)
    gwas_df.loc[gwas_df['sample_size'].isnull(), 'sample_size'] = median_sample_size
  else:
    gwas_df.loc[gwas_df['sample_size'].isnull(), 'sample_size'] = sample_size_d[phenotype]['sample_size']

###############################################################


################################
##   FILLING WITH FREQUENCY   ##
################################

# RUN ONLY ONCE, NOT SURE IF WORKS IN ITS CURRENT FORM
def build_maf_dictionary(genotype_folder, genotype_metadata, pkl_file):
  data = Parquet.study_from_parquet(genotype, genotype_metadata).variant_metadata[["id", "allele_1_frequency"]]
  variants_id = data["id"]
  mafs = data["allele_1_frequency"]

  import pickle
  with open(pkl_file, "wb") as dd:
    pickle.dump(dict(zip(list(variants_id), list(mafs))), dd)


# WORKS
regex = re.compile("chr.*_.*_(.*)_(.*)_b38")
def extract_alleles(panel_variant_id):
  gtex_allele_1, gtex_allele_2 = regex.search(panel_variant_id).groups()
  return gtex_allele_1, gtex_allele_2


def is_maf_available(comps):
  return not numpy.isnan(comps["frequency"])


def fill_with_gtex_freq(comps, snp_dictionary):

  variant_id = comps["panel_variant_id"]
  try:
    gtex_allele_1, gtex_allele_2 = extract_alleles(variant_id)
  except:
    pass
 
  if not is_maf_available(comps):
    try:
      maf = snp_dictionary[variant_id]
      if comps['effect_allele'] == gtex_allele_1 and comps['non_effect_allele'] == gtex_allele_2:
        comps['frequency'] = 1 - maf
      elif comps['non_effect_allele'] == gtex_allele_1 and comps['effect_allele'] == gtex_allele_2:
        comps['frequency'] = maf
    except:
      comps['frequency'] = "NA"
  return comps


def load_pkl(pkl_file):
  import pickle
  with open(pkl_file, "rb") as pkl:
    dd = pickle.load(pkl)
  return(dd)


def gtex_freq(comps, snp_dictionary):
  variant_id = comps["panel_variant_id"]
  try:
    gtex_allele_1, gtex_allele_2 = extract_alleles(variant_id)
  except:
    pass

  if not is_maf_available(comps):
    try:
      maf = snp_dictionary[variant_id]
      if comps['effect_allele'] == gtex_allele_1 and comps['non_effect_allele'] == gtex_allele_2:
        return(1-maf)
      elif comps['non_effect_allele'] == gtex_allele_1 and comps['effect_allele'] == gtex_allele_2:
        return(maf)
    except:
      return("NA")

  return(comps['frequency'])


### PATCH GWAS ###
def patch_gwas(gwas_df, mafs, phenotype, sample_sizes):
  def _fill_with_gtex_freq(row):
    return fill_with_gtex_freq(row, mafs)

  mafs_lst = []
 
  logging.info("Filling with missing MAF's.")
 
  for i, line in gwas_df.iterrows():
    mafs_lst.append(gtex_freq(line, mafs))
    if i % 100000 == 0: 
      print(i)
  
  gwas_df.loc[:,'frequency'] = mafs_lst
  
  logging.info("MAF's filled in.")
  logging.info("Filling with sample size.")

  # embed()
  if sample_sizes[phenotype]['binary']:
    gwas_df['n_cases'] = sample_sizes[phenotype]['cases']
  else:
    gwas_df['n_cases'] = "NA"

  complete_missing_sample_sizes(gwas_df, phenotype, sample_sizes)
  
  logging.info("Sample sizes filled in.")
  return(gwas_df)

##################################

def create_job_file(pheno_name, gwas_filename, jobs_folder, logs_folder, resource_manager, HOME=None):

  if HOME is None:
    if resource_manager.lower() == "slurm":
      HOME = "/home/rodrigo"
    elif resource_manager.lower() == "pbs" or resource_manager.lower() == "torque":
      HOME = "/group/im-lab/nas40t2/rbonazzola" 

  print("home directory: ", HOME)
 
  def write_job_content(job_file, resource_manager):
    ## SLURM ## Intended mainly for GCP
    if resource_manager.lower() == "slurm":
      job_file.write("""#!/bin/bash\n\
#SBATCH --job-name={job_name}\n\
#SBATCH --time=5:00:00\n\
#SBATCH --mem=30gb\n\
#SBATCH --cpus-per-task=8\n\
#SBATCH --output={log}\n\
#SBATCH --error={err}\n\n\

if [ -z $SLURM_SUBMIT_DIR ]; then cd $SLURM_SUBMIT_DIR; fi

source activate py3
python /home/rodrigo/software/genomic_tools/src/genomic_tools/patch_imputed.py --gwas_input_filename {gwas}

""".format(job_name=pheno_name + "_patch_imputed", \
           gwas=gwas_filename, \
           log=os.path.join(logs_folder, pheno_name + ".log"), \
           err=os.path.join(logs_folder, pheno_name + ".err")))
    ## TORQUE ## Intended mainly for CRI
    elif resource_manager.lower() == "pbs" or resource_manager.lower() == "torque":
      job_file.write("""#!/bin/bash\n\
#PBS -N {job_name}\n\
#PBS -l walltime=5:00:00\n\
#PBS -l mem=30GB\n\
#PBS -o {log}\n\
#PBS -e {err}\n\n\

if [ -z $PBS_O_WORKDIR ]; then cd $PBS_O_WORKDIR; fi

module load gcc/6.2.0
module load python/3.5.3

python3 {execut} --gwas_input_filename {gwas}
""".format(execut=os.path.join(HOME, "GTEx/v8/summary-imputation/patching/patch_imputed.py"),
           job_name=pheno_name + "_patch_imputed", \
           gwas=gwas_filename, \
           log=os.path.join(logs_folder, pheno_name + ".log"), \
           err=os.path.join(logs_folder, pheno_name + ".err")))


  if not os.path.exists(jobs_folder):
    os.makedirs(jobs_folder)

  if not os.path.exists(logs_folder):
    os.makedirs(logs_folder)

  job_filename = os.path.join(jobs_folder, pheno_name + "_patch_imputed.sh")

  with open(job_filename, "w") as job_file:
    write_job_content(job_file, resource_manager)

  return(job_filename)

################################################################


def get_phenotype_from_gwas_name(gwas_filename):
  phenoname = gwas_filename.split(".txt")[0].split("imputed_")[-1]
  return(phenoname)


def submit_job(job, resource_manager):
  if resource_manager.lower() == "slurm":
    subprocess.call(["sbatch",job])
  elif resource_manager.lower() == "pbs" or resource_manager.lower() == "torque":
    subprocess.call(["qsub", job])
  else:
    logging.error("Use a valid resource manager (Slurm or PBS/Torque.)")
    stop()


def run(args):
  
  logging.info("Reading sample size...")
  sample_sizes = read_sample_size(args.sample_size_file)

  if not os.path.exists(args.imputed_results_folder):
    os.makedirs(args.imputed_results_folder)

  if args.create_jobs_for_cluster:
  
    if args.fetch_gwas_from_gcs:
      gwasfile_lst = [os.path.split(x)[1] for x in subprocess.check_output(shlex.split("gsutil ls gs://summary_stats_imputation/processed_summary_imputation_deprecated")).decode().strip().split() if "txt" in x] 
    else:
      gwasfile_lst = os.listdir(args.imputed_results_folder)

    if not os.path.exists(args.jobs_folder):
      os.makedirs(args.jobs_folder)
    if not os.path.exists(args.logs_folder):
      os.makedirs(args.logs_folder)

    for gwas_filename in gwasfile_lst:
      pheno_name = get_phenotype_from_gwas_name(gwas_filename)
      # gwas_filename = os.path.join(args.imputed_results_folder, gwas_filename)
      if pheno_name in sample_sizes.keys():
        logging.info("Creating job for: {pheno}".format(pheno=pheno_name))
        job_filename = create_job_file(pheno_name, gwas_filename, args.jobs_folder, args.logs_folder, args.resource_manager)
        if args.submit_jobs_to_cluster:
          submit_job(job_filename, args.resource_manager)

  else:
    # the "then" part calls the "else" part
    if not os.path.exists(args.imputed_results_folder):
      os.makedirs(args.imputed_results_folder)

    imputed_gwas_path = os.path.join(args.imputed_results_folder, args.gwas_input_filename)

    if args.fetch_gwas_from_gcs:
      if not os.path.exists(imputed_gwas_path):
        print("Copying data from Google Cloud Storage...")
        subprocess.call(shlex.split("gsutil cp gs://summary_stats_imputation/processed_summary_imputation_deprecated/{gwas} {imp_gwas_path}".format(imp_gwas_path=imputed_gwas_path, gwas=args.gwas_input_filename)))
      else:
        print("File already downloaded in the past.")
  
    logging.info("Loading MAF dictionary...")
    mafs = load_pkl(args.maf_dictionary)
    
   
    """
    print("Filtering for MAF > 0.01...")
    mafs = {snp:mafs[snp] for snp in mafs.keys() if mafs[snp] > 0.01}
    
    with open("gtex_variants_maf_filtered_0.01.pkl", "wb") as pkl:
      pickle.dump(mafs, pkl) 
    print(len(mafs))
    """

    pheno_name = get_phenotype_from_gwas_name(args.gwas_input_filename)
    logging.info("Loading GWAS for phenotype {}".format(pheno_name))
    gwas_df = read_imputed_gwas(os.path.join(args.imputed_results_folder,args.gwas_input_filename))

    # main function
    logging.info("Trying to patch GWAS...")
    gwas_df = patch_gwas(gwas_df, mafs, pheno_name, sample_sizes)

    if args.gwas_output_filename is None:
      gwas_output_filename = os.path.join(args.gwas_output_folder, os.path.split(args.gwas_input_filename)[1])
    else:
      gwas_output_filename = os.path.join(args.gwas_output_folder, args.gwas_output_filename)

    logging.info(gwas_output_filename)
    logging.info("Saving patched GWAS in {}...".format(gwas_output_filename))

    if not os.path.exists(args.gwas_output_folder):
      os.makedirs(args.gwas_output_folder)

    if gwas_output_filename.endswith(".gz"):
      gwas_df.to_csv(gwas_output_filename, sep="\t", compression = "gzip", na_rep = "NA", index = False)
    else:
      gwas_df.to_csv(gwas_output_filename, sep="\t", index = False, na_rep = "NA")

    if args.save_in_gcs:
      logging.info("Copying data to GCS bucket")
      subprocess.call(shlex.split("gsutil cp {} gs://summary_stats_imputation/processed_summary_imputation".format(gwas_output_filename)))

    logging.info("Done.")


if __name__ == "__main__":

  import argparse
  
  parser = argparse.ArgumentParser()

  parser.add_argument("--fetch_gwas_from_gcs", action="store_true", default = False)
  parser.add_argument("--imputed_results_folder", default = "/group/im-lab/nas40t2/rbonazzola/GTEx/v8/summary-imputation/collection/processed_summary_imputation/")
  parser.add_argument("--sample_size_file", default="/group/im-lab/nas40t2/rbonazzola/GTEx/v8/summary-imputation/patching/sample_sizes.txt")
  parser.add_argument("--sample_size_available", action="store_true", default = False)
  parser.add_argument("--maf_dictionary", default="/group/im-lab/nas40t2/rbonazzola/GTEx/v8/summary-imputation/patching/gtex_variant_mafs.pkl")

  parser.add_argument("--create_jobs_for_cluster", action="store_true", default=False)
  parser.add_argument("--submit_jobs_to_cluster", action="store_true", default=False)
  parser.add_argument("--resource_manager", default="SLURM")
  parser.add_argument("--jobs_folder", default="/home/rodrigo/imputation/patched_imputation/jobs")
  parser.add_argument("--logs_folder", default = "/home/rodrigo/imputation/patched_imputation/logs")

  parser.add_argument("--gwas_input_filename", default="/mnt/disks/data/gwas/UKB_50_Standing_height.txt.gz") # TO TEST
  parser.add_argument("--gwas_output_folder", default="/group/im-lab/nas40t2/rbonazzola/GTEx/v8/summary-imputation/patching/results")
  parser.add_argument("--gwas_output_filename", default=None)
  parser.add_argument("--save_in_gcs", default=None)

  parser.add_argument("--genotype_folder", default = "/mnt/disks/data/gtex_genotype")
  parser.add_argument("--genotype_metadata", default = "gtex_v8_eur_itm.variants_metadata.parquet")

  args = parser.parse_args()

  run(args)
