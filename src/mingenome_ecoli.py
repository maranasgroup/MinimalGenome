from codebase.core import *
from codebase.util import *
import pandas as pd

project_name = 'Ecoli'
GSM_name = "/ijo1366.xml"
lpfile_name = "mingenome_Ecoli.lp"

data_dir = "./data/" + project_name
GSM_dir = "./data/" + project_name + GSM_name
sequence_dir = "./data/" + project_name + "/NC_000913.3"
offrxn_f = None

### setup input #####
# util = utility(sequence_dir)
# util.getAllGenesfromGB()
# util.get_promoter_end_loc()
# util.move_start_site_to_back(data_dir)
# util.exchange_promoter_key_value(data_dir)
# util.map_genes_between_models(data_dir)

# ### run MinGenome algorithm #####
mg = MinGenome(GSM_dir)
model = mg.modelCobra(data_dir,offrxn_f)

eg_f = "./data/" + project_name + "/essentialGene.txt"
# eg_f = "./data/" + project_name + "/essentialGenes_Daniel.txt"
parameters_f = "./data/" + project_name + "/genes_and_promoters.xlsx"
# reg_f = "./data/" + project_name + "/regulatorGene.txt"
reg_f = None
TU_Json_file="./data/" + project_name + "/gene_promoter_dict.json"
lpfilename="./out/" + lpfile_name

mg.build_MIP_by_Cobrapy(model, 0.98,
		eg_f,parameters_f,reg_f,TU_Json_file,lpfilename)