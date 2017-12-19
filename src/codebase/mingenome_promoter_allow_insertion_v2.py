import pandas
import pdb
import re
import json
import cobra
import pulp
import pdb
import cPickle as pickle
import itertools
# from framed.cobra.simulation import FBA

# --------- analysis result
def getAllgenesInResult():
    LD = "/Users/linuswang/Box Sync/Research@PennState/My Paper/LD.xlsx"
    all = "../data/e_coli/genes_and_promoters.xlsx"
    df_all = pandas.read_excel(all,sheetname='all_clear_v2',\
            index_col='gene_or_promoter')
    df_LD = pandas.read_excel(LD,sheetname='result_with_promoters',\
            index_col='rank')
    dict_LD = df_LD.to_dict(orient='index')

    result = pandas.DataFrame()
    for rank, value in dict_LD.iteritems():
        start = value['start_gene']
        end = value['end_gene']
        sub_df_all = df_all[start:end]
        # pdb.set_trace()
        result = result.append(sub_df_all)
    result = result[result['class']=='gene']
    result.to_csv('LD_all_genes.csv')
    pdb.set_trace()    

def getOverlap(a, b):
     return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def get_overlap():

    all = "../data/e_coli/genes_and_promoters.xlsx"
    df = pandas.read_excel(all,sheetname='all')
    # pdb.set_trace()
    for i in range(0, len(df)):
        a = df.iloc[i]['start'], df.iloc[i]['end']
        for j in range(i+1, i+7):
            b = df.iloc[j]['start'], df.iloc[j]['end']
            if getOverlap(a, b) != 0:
                print df.iloc[i]['gene_or_promoter'], df.iloc[j]['gene_or_promoter']

def find_genes_in_TU():
    f = "../data/e_coli/3556_Transcription-Units_from_a_list (3).txt"
    pattern = re.compile('b[0-9]+')
    with open(f,'rb') as rf:
        lines = rf.readlines()
        for line in lines:
            print pattern.findall(line)

# some genes or TUs have more than one promoters, keep all of them or delete all
# at the same time. As a result, 
def remove_duplicate_promoter():
    all = "../data/e_coli/genes_and_promoters.xlsx"
    df = pandas.read_excel(all,sheetname='all_clear')
    df_promoters = df[df['class']=='promoter']

    # pdb.set_trace()
    for i in range(0, len(df)):
        to_delete = {}
        a_name = df.iloc[i]['gene_or_promoter']
        a_genes = df.iloc[i]['genes_in_TU']
        a_start = df.iloc[i]['start']
        a_end = df.iloc[i]['end']
        to_delete[a_name] = [a_start, a_end]

        # if a_name == 'PM00014': pdb.set_trace()
        for j in range(i+1, len(df)):
            b_genes = df.iloc[j]['genes_in_TU']
            if a_genes != b_genes: break
            b_name = df.iloc[j]['gene_or_promoter']
            b_start = df.iloc[j]['start']
            b_end = df.iloc[j]['end']
            to_delete[b_name] = [b_start, b_end]

        keep_promoter = a_name
        for gene,loc in to_delete.iteritems():
            if loc[0] <= a_start and loc[1] >= a_end:
                a_start = loc[0]
                a_end = loc[1]
                keep_promoter = gene
        for gene in to_delete.keys():
            if gene != keep_promoter:
                print gene
def get_promoter_removed_df():
    duplicate_promoter = pandas.read_csv('./promoter_to_delete.txt',\
                            header=None)[0].tolist()
    # pdb.set_trace()
    all = "../data/e_coli/genes_and_promoters.xlsx"
    df = pandas.read_excel(all,sheetname='all_clear')

    non_dup_df = df[~df['gene_or_promoter'].isin(duplicate_promoter)]

    non_dup_df.to_csv('../data/e_coli/genes_and_promoters_v2.csv')

def move_start_site_to_back():
    all = "../data/e_coli/genes_and_promoters.xlsx"
    df = pandas.read_excel(all,sheetname='all_clear_v2')

    start_list = []
    for i in range(5, len(df)):
        to_delete = {}
        a_name = df.iloc[i]['gene_or_promoter']
        a_genes = df.iloc[i]['genes_in_TU']
        a_start = df.iloc[i]['start']
        a_end = df.iloc[i]['end']

        # if a_name == 'PM0-8866': pdb.set_trace()
        # if a_name == 'PM00014': pdb.set_trace()
        for j in range(i-6, i):
            b_end = df.iloc[j]['end']
            if b_end > a_start:
                a_start = b_end
        start_list.append(a_start)
    pandas.DataFrame(start_list).to_csv("../data/e_coli/genes_and_promoters_new_start.csv")

def exchange_promoter_key_value():
    all = "../data/e_coli/genes_and_promoters.xlsx"
    df = pandas.read_excel(all,sheetname='all_clear')
    # df_promoters = df[df['class']=='promoter'][['gene_or_promoter','genes_in_TU']]

    # get the dictionary of promoters -> list of genes
    promoters = df[df['class']=='promoter'][['gene_or_promoter']]['gene_or_promoter'].tolist()
    genes_str = df[df['class']=='promoter'][['genes_in_TU']]['genes_in_TU'].tolist()
    new_genes = []
    for gene in genes_str:
        new_genes.append(gene.replace('[','').replace(']','').split())

    promoter_genes_dict = {}
    for idx,promoter in enumerate(promoters):
        promoter_genes_dict[promoter] = new_genes[idx]

    # change to gene to list of promoters
    gene_promoter_dict = {}
    for idx,genes in enumerate(new_genes):
        for gene in genes:
            if gene not in gene_promoter_dict.keys():
                gene_promoter_dict[gene] = [promoters[idx]]
            else:
                gene_promoter_dict[gene].append(promoters[idx])
    with open('../data/e_coli/gene_promoter_dict.json', 'w') as fp:
        json.dump(gene_promoter_dict, fp)


def modelCobra(data_dir):
    model = cobra.io.read_sbml_model(data_dir)
    # pdb.set_trace()
    offrxn_f = "../data/e_coli/off_rxn_glucose_aerobic.txt"
    # eg_f = "../data/e_coli/essentialGene_with_Buchenera.txt"
    offrxn_df = pandas.read_csv(offrxn_f,header=None)
    # essential_genes[0] = map(lambda x: x.upper(), essential_genes[0])

    # model.optimize()
    # print model.solution.status
    # print model.solution.f
    # model.summary()
    # off_rxn = offrxn_df[0].tolist()
    # for r in off_rxn:
    #     r = model.reactions.get_by_id(r)
    #     r.lower_bound = 0
    #     r.upper_bound = 0
    # for r in model.reactions:
    #     print r.id,r.lower_bound,r.upper_bound

    model.optimize()
    print model.solution.status
    print model.solution.f
    model.summary()
    return model

def get_S(model):
    """build the stoichiometric matrix at a specific growth rate"""
        # intialize to 0
    # S = dok_matrix((len(self.metabolites), len(self.reactions)))
    S = {}
    # populate with stoichiometry
    for i, r in enumerate(model.reactions):
        for met, value in r._metabolites.iteritems():
            #met_index = self.metabolites.index(met)
            if met.id not in S:
                S[met.id] = {}
            if hasattr(value, "subs"):
                S[met.id][r.id] = float(value.subs(mu, growth_rate))
            else:
                S[met.id][r.id] = float(value)
    return S

def build_MIP_by_Cobrapy(me,mu):

    M = 1000
    
    ##### sets
    # TU
    with open('../data/e_coli/gene_promoter_dict.json') as data_file:    
        TUs = json.load(data_file)
    # essential genes
    # eg_f = "../data/e_coli/essentialGene.txt"
    eg_f = "../data/e_coli/essentialGene_with_Buchenera.txt"
    essential_genes = pandas.read_csv(eg_f,header=None)
    # essential_genes[0] = map(lambda x: x.upper(), essential_genes[0])
    essential_genes[0] = "u_G_" + essential_genes[0].astype(str)
    essential_genes = essential_genes[0].tolist()

    # regulator genes
    reg_f = "../data/e_coli/regulatorGene.txt"
    reg_genes = pandas.read_csv(reg_f,header=None)
    # essential_genes[0] = map(lambda x: x.upper(), essential_genes[0])
    reg_genes[0] = "u_G_" + reg_genes[0].astype(str)
    reg_genes = reg_genes[0].tolist()

    ### parameters
    # f = "../data/e_coli/allgenes_e_coli.txt"
    f = "../data/e_coli/genes_and_promoters.xlsx"
    df = pandas.read_excel(f,sheetname='all_clear_v2')

    test_all_genes = df["gene_or_promoter"].tolist()
    not_shared = []
    for gene in TUs.keys():
        if gene not in test_all_genes:
            not_shared.append(gene)
    # print [gene for gene in TUs.keys() and gene not in test_all_genes]
    # print not_shared
    # assert(set(TUs.keys()).issubset(test_all_genes))

    df["gene_or_promoter"] = "u_G_" + df["gene_or_promoter"].astype(str)
    no_start = df[df['cannot_as_start']==1]["gene_or_promoter"].tolist()

    genes = df["gene_or_promoter"].tolist()
    
    end = df[['gene_or_promoter','start']].set_index('gene_or_promoter')\
                    .T.to_dict('list')
    start = df[['gene_or_promoter','start_if_select_as_start']]\
                    .set_index('gene_or_promoter').T.to_dict('list')
    reactions = [r_id.id for r_id in me.reactions]
    metabolites = [m_id.id for m_id in me.metabolites]
    
    #### variables
    v = pulp.LpVariable.dicts("v", reactions, 0, M, cat='Continuous')
    x = pulp.LpVariable.dicts("x", genes, cat='Binary')
    y = pulp.LpVariable.dicts("y", genes, cat='Binary')
    z = pulp.LpVariable.dicts("z", genes, cat='Binary') # can be defined as continuous
    w = pulp.LpVariable.dicts("w", genes, cat='Binary') # can be defined as continuous

    # print len(x), len(y),len(z),len(v)
    #### define model
    lp_prob = pulp.LpProblem("MaxDeletion", pulp.LpMaximize)
    
    #### objective
    lp_prob += (pulp.lpSum([y[j]*end[j][0] for j in genes]) 
              - pulp.lpSum([x[j]*start[j][0] for j in genes])), "Max_length"
    
                    
    def addReactionIndicator(lp_prob):
        for r in me.reactions:
            rgenes = r.genes
            GPR = r.gene_reaction_rule
            if 's0001' in GPR: continue

            # pdb.set_trace()
            # no genes
            if len(rgenes) == 0:
                continue
            # single gene
            # if ('and' and 'AND' and 'or' and 'OR') not in GPR:
            if 'and' not in GPR \
            and 'AND' not in GPR \
            and 'or' not in GPR \
            and 'OR' not in GPR:
                # print GPR, genes
                assert(len(rgenes) == 1)
                for gene in rgenes:
                    gene = gene.id
                    label = "knockout" + str(gene)
                    gene = "u_G_" + gene
                    lp_prob += v[r.id] - (1-z[gene])*M <= 0, label + "_UB_" + r.id
                    lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, label + "_LB_" + r.id

            # enzyme complex
            elif (('and' or 'AND') in GPR) and  (('or' or 'OR') not in GPR):
                # print genes
                # print GPR
                assert(len(rgenes) > 1)
                for gene in rgenes:
                    gene = gene.id
                    label = "knockout_" + str(gene)
                    gene = "u_G_" + gene
                    lp_prob += v[r.id] - (1-z[gene])*M <= 0, label + "_UB_" + r.id
                    lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, label + "_LB" + r.id

            # isozymes
            elif (('and' or 'AND') not in GPR) and  (('or' or 'OR') in GPR):
                # print GPR
                lp_prob += v[r.id] - M <= 0, "knockout" + r.id + "Ori_UB"
                lp_prob += v[r.id] - (-M) >= 0, "knockout" + r.id + "Ori_LB"
                assert(len(rgenes) > 1)
                lp_prob += v[r.id] - M * pulp.lpSum(1-z['u_G_'+j.id] for j in rgenes) <=0, "knockout" + r.id + '_UB'
                lp_prob += v[r.id] - (-M) * pulp.lpSum(1-z['u_G_'+j.id] for j in rgenes) >=0, "knockout" + r.id + '_LB'
            # more complicated GPRs
            else:
                # print r.id
                # print GPR
                proteins = [protein.replace("( ","").replace(" )","").split(' and ') \
                                     for protein in GPR.split(' or ')]
                
                commonGenes = set(proteins[0])
                # if len(gpr.proteins) > 1:
                for protein in proteins[1:]:
                    commonGenes.intersection_update(protein)
                nonCommonGenesList = []
                for protein in proteins:
                    nonCommonGenes = []
                    for gene in protein:
                        if gene not in commonGenes:
                            nonCommonGenes.append(gene)
                    nonCommonGenesList.append(nonCommonGenes)
                    
                for gene in commonGenes:
                    # gene = gene.id
                    label = "knockout" + str(gene)
                    gene = "u_G_" + gene
                    lp_prob += v[r.id] - (1-z[gene])*M <= 0, label + "_UB_" + r.id
                    lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, label + "_LB_" + r.id
                    
                allCombination = list(itertools.product(*nonCommonGenesList))
                # print allCombination
                for i,genesC in enumerate(allCombination):
                    lp_prob += v[r.id] - M * pulp.lpSum(1-z['u_G_'+j] for j in genesC) <=0,\
                        "knockout" + r.id + '_UB_' + str(i)
                    lp_prob += v[r.id] - (-M) * pulp.lpSum(1-z['u_G_'+j] for j in genesC) >=0,\
                        "knockout" + r.id + '_LB_' + str(i)

    #### constraint
    print "add reaction indicator"    
    addReactionIndicator(lp_prob)

    # ME model constraint
    S = get_S(me) # growth rate is 0.3
    
    print "add GSM constraint"    
    for i in metabolites:
        label = "mass_balance_%s"%i
        dot_S_v = pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()])
        condition = dot_S_v == 0
        lp_prob += condition, label

    print "add cutting constraints"    
    # cut in the genome        
    lp_prob += pulp.lpSum(y[j] for j in genes) == 1, "end"
    lp_prob += pulp.lpSum(x[j] for j in genes) == 1, "start"

    # cut genes between start and end
    # for i,gene in enumerate(genes):
    #     lp_prob += pulp.lpSum(x[j] for j in genes[0:i+1]) - pulp.lpSum(y[j] for j in genes[0:i+1]) - z[gene] == 0,'indicator' + str(gene)
    
    A = pulp.LpAffineExpression()

    for i,gene in enumerate(genes):
        # print i
        A.addterm(x[gene],1) #pulp.lpSum(x[j] for j in genes[0:i+1])
        A.addterm(y[gene],-1)
        # lp_prob += A - z[gene] == 0,'indicator' + str(gene)
        lp_prob += A - w[gene] == 0,'indicator' + str(gene)
        lp_prob += w[gene] >= z[gene],'indictart_1_'+str(gene)
    lp_prob += pulp.lpSum(w[j] for j in genes) - pulp.lpSum(z[j] for j in genes) == 1, "allow_one_insertion"


    # TUs 
    print "add TU constraint"
    for gene,promoters in TUs.iteritems():
        if gene in not_shared: continue
        gene = 'u_G_' + gene
        len_pro = len(promoters)
        lp_prob += z[gene] - pulp.lpSum(z['u_G_'+j] for j in promoters) + \
                    (len_pro - 1) >= 0,'TU_all_'+gene
        for pro in promoters:
            pro = 'u_G_' + pro
            lp_prob += z[gene] - z[pro] <=0, 'TU_'+gene+'_'+pro

    # some overlapped region cannot be selected as the start of deletions
    print "add no start and essential genes"
    for gene in no_start:
        lp_prob += x[gene] == 0, 'no_start_'+gene

    # # knock out transcription of cutted genes
    # for gene in genes:
    #     label = "knockout" + str(gene)
    #     # pdb.set_trace()
    #     if gene in v.keys():
    #         lp_prob += v[gene] - (1-z[gene])*M <= 0, label
    
    # add essential genes that cannot be deleted
    for eg in essential_genes:
        if eg in genes:
            lp_prob += z[eg] == 0
    # add regulation genes that cannot be deleted
    for eg in reg_genes:
        # remove the part joint with essential genes
        if (eg in genes) and (eg not in essential_genes):
            lp_prob += z[eg] == 0

    # reaction bounds
    for r in me.reactions:
        # (lb,up) = me.bounds[r_id]
        v[r.id].lowBound = r.lower_bound
        v[r.id].upBound = r.upper_bound

    # specify growth rate
    # lp_prob += v['biomass_dilution'] == mu, "growth_rate"
    # v['biomass_dilution'].lowBound = mu
    # v['biomass_dilution'].upBound = mu
    v['BIOMASS_Ec_iJO1366_WT_53p95M'].lowBound = mu
    # v['BIOMASS_Ec_iJO1366_WT_53p95M'].upBound = mu

    v['BIOMASS_Ec_iJO1366_core_53p95M'].lowBound = 0
    v['BIOMASS_Ec_iJO1366_core_53p95M'].upBound = 0

    # lp file is too larget to write
    # lp_prob.writeLP("mingenome_ecoli_with_regulation.lp")
    # lp_prob.writeLP("mingenome_ecoli_with_regulation_Bun.lp")
    lp_prob.writeLP("mingenome_ecoli_with_regulation_Bun_allowinsertion1.lp")

    # options = [epgap, epagap, epint, epopt, eprhs]
    
    # pulp_solver = pulp.solvers.CPLEX_CMD(path=None, keepFiles=0, mip=1, msg=1, \
    #     options=['mip tolerances mipgap 0', 'mip tolerances absmipgap 0', 'mip tolerances integrality 0', 'simplex tolerances optimality 1E-9','simplex tolerances feasibility 1E-9',], timelimit=1200)
    # lp_prob.solve(pulp_solver)
    # print("Status:", pulp.LpStatus[lp_prob.status])

    # pdb.set_trace()
    # with open("./maxlength_solution",'wb') as f:
    #     f.write("startus: " + str(pulp.LpStatus[lp_prob.status]) + '\n')
    #     # f.write("Objctive_max_deletion = " + value(lp_prob.objective) + '\n')
    #     for v in lp_prob.variables():
    #         f.write(str(v.name) + "=" + str(v.varValue) + '\n')

if __name__ == '__main__':
    # find_genes_in_TU()
    # get_overlap()
    # remove_duplicate_promoter()
    # get_promoter_removed_df()
    # move_start_site_to_back()
    # exchange_promoter_key_value()
    
    # #### solve
    data_dir = "../data/e_coli/iJO1366.xml"
    model = modelCobra(data_dir)
    mip = build_MIP_by_Cobrapy(model, 0.98)

    #### analyze result
    # getAllgenesInResult()
