import pandas
import pdb
import re
import json
import cobra
import pulp
import pdb
import cPickle as pickle
import itertools

class MinGenome(object):
    """docstring for MinGenome"""

    def __init__(self, arg):
        super(MinGenome, self).__init__()
        self.model_file = arg
    
    def modelCobra(self,data_dir,
        offrxn_f=''):

        model = cobra.io.read_sbml_model(self.model_file)
        # model = cobra.io.load_matlab_model(self.model_file)
        if offrxn_f != None:
            offrxn_df = pandas.read_csv(offrxn_f,header=None)
            # add bound to constrain the flux of these reactions to be zero

        model.optimize()
        print model.solution.status
        print model.solution.f
        model.summary()

        # df_medium = pandas.read_excel(data_dir+'/medium.xls',sheetname="LB medium")
        # exchange_rxns = df_medium['Reaction No.'].tolist()
        # LBs = df_medium['LB'].tolist()
        # UBs = df_medium['UB'].tolist()

        # for idx,r in enumerate(exchange_rxns):
        #     rxn = model.reactions.get_by_id(r)
        #     rxn.lower_bound = LBs[idx]
        #     rxn.upper_bound = UBs[idx]

        # model.optimize()
        # print model.solution.status
        # print model.solution.f
        # model.summary()
        return model
    
    def modelCobra_change_GPR(self,data_dir,oldmodel_dir,
        offrxn_f=''):

        oldmodel = cobra.io.read_sbml_model(oldmodel_dir)
        model = cobra.io.load_matlab_model(self.model_file)
        # model = cobra.io.read_sbml_model(self.model_file)
        for r in model.reactions:
            if "biomass" in r.id: continue
            oldGPR = oldmodel.reactions.get_by_id(r.id).gene_reaction_rule
            r.gene_reaction_rule = oldGPR
        if offrxn_f != None:
            offrxn_df = pandas.read_csv(offrxn_f,header=None)
            # add bound to constrain the flux of these reactions to be zero

        model.optimize()
        print model.solution.status
        print model.solution.f
        model.summary()

        df_medium = pandas.read_excel(data_dir+'/medium.xls',sheetname="LB medium")
        exchange_rxns = df_medium['Reaction No.'].tolist()
        LBs = df_medium['LB'].tolist()
        UBs = df_medium['UB'].tolist()

        for idx,r in enumerate(exchange_rxns):
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound = LBs[idx]
            rxn.upper_bound = UBs[idx]

        model.optimize()
        print model.solution.status
        print model.solution.f
        model.summary()
        return model

    def get_S(self,model):
        """build the stoichiometric matrix at a specific growth rate"""
        S = {}
        # populate with stoichiometry
        for i, r in enumerate(model.reactions):
            # if r.id == 'R00026': pdb.set_trace()
            for met, value in r._metabolites.iteritems():
                #met_index = self.metabolites.index(met)
                # if met.id == "C00185_c": pdb.set_trace()
                if met.id not in S:
                    S[met.id] = {}
                    S[met.id][r.id] = float(value)
                else:
                    S[met.id][r.id] = float(value)
        return S

    def build_MIP_by_Cobrapy(self,me,mu,
        eg_f = "./data/e_coli/essentialGene.txt",
        parameters_f = "./data/e_coli/genes_and_promoters.xlsx",
        reg_f = "./data/e_coli/regulatorGene.txt",
        TU_Json_file='./data/e_coli/gene_promoter_dict.json',
        lpfilename="./out/mingenome_ecoli_with_regulation_Bun.lp"):

        M = 1000
        
        ############# sets ################################
        # TU
        with open(TU_Json_file) as data_file:    
            TUs = json.load(data_file)
        
        # essential genes        
        essential_genes = pandas.read_csv(eg_f,header=None)
        essential_genes[0] = "u_G_" + essential_genes[0].astype(str)
        essential_genes = essential_genes[0].tolist()

        # regulator genes
        if reg_f != None:
            reg_genes = pandas.read_csv(reg_f,header=None)
            reg_genes[0] = "u_G_" + reg_genes[0].astype(str)
            reg_genes = reg_genes[0].tolist()

        ############# parameters ################################       
        df = pandas.read_excel(parameters_f,sheetname='all_clear_v2')

        test_all_genes = df["gene_or_promoter"].tolist()
        not_shared = []
        for gene in TUs.keys():
            if gene not in test_all_genes:
                not_shared.append(gene)

        df["gene_or_promoter"] = "u_G_" + df["gene_or_promoter"].astype(str)
        no_start = df[df['cannot_as_start']==1]["gene_or_promoter"].tolist()

        genes = df["gene_or_promoter"].tolist()
        
        end = df[['gene_or_promoter','start']].set_index('gene_or_promoter')\
                        .T.to_dict('list')
        start = df[['gene_or_promoter','start_if_select_as_start']]\
                        .set_index('gene_or_promoter').T.to_dict('list')

        reactions = [r_id.id for r_id in me.reactions]
        metabolites = [m_id.id for m_id in me.metabolites]
        
        ############# variables ################################
        v = pulp.LpVariable.dicts("v", reactions, 0, M, cat='Continuous')
        x = pulp.LpVariable.dicts("x", genes, cat='Binary')
        y = pulp.LpVariable.dicts("y", genes, cat='Binary')
        z = pulp.LpVariable.dicts("z", genes, cat='Binary') 
        # z can be defined as continuous


        ############# define model ################################
        lp_prob = pulp.LpProblem("MaxDeletion", pulp.LpMaximize)
        
        ############# objective ################################
        lp_prob += (pulp.lpSum([y[j]*end[j][0] for j in genes]) 
                  - pulp.lpSum([x[j]*start[j][0] for j in genes])), "Max_length"
        
                        
        def addReactionIndicator(lp_prob):
            for r in me.reactions:
                rgenes = r.genes
                GPR = r.gene_reaction_rule
                GPR = GPR.replace('\n','')
                GPR = GPR.replace('__10__','')
                # GPR = GPR.replace('Bsu1823a','BSU18239')
                # GPR = GPR.replace('Bsu3567','BSU35699')
                # print rgenes
                # print GPR

                if 's0001' in GPR: continue # not mapped gene in iJO1366
                if 'BG12900' in GPR: continue # not mapped gene in iYO844

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
                        gene = gene.id.replace('__10__','')
                        label = "knockout" + str(gene)
                        gene = "u_G_" + gene
                        lp_prob += v[r.id] - (1-z[gene])*M <= 0, \
                                            label + "_UB_" + r.id
                        lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, \
                                            label + "_LB_" + r.id

                # enzyme complex
                elif (('and' or 'AND') in GPR) and  (('or' or 'OR') not in GPR):
                    # print genes
                    # print GPR
                    assert(len(rgenes) > 1)
                    for gene in rgenes:
                        gene = gene.id.replace('__10__','')
                        label = "knockout_" + str(gene)
                        gene = "u_G_" + gene
                        lp_prob += v[r.id] - (1-z[gene])*M <= 0, \
                                            label + "_UB_" + r.id
                        lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, \
                                            label + "_LB" + r.id

                # isozymes
                elif (('and' or 'AND') not in GPR) and  (('or' or 'OR') in GPR):
                    # print GPR
                    lp_prob += v[r.id] - M <= 0, "knockout" + r.id + "Ori_UB"
                    lp_prob += v[r.id] - (-M) >= 0, "knockout" + r.id + "Ori_LB"
                    assert(len(rgenes) > 1)
                    lp_prob += v[r.id] - M * pulp.lpSum(1-z['u_G_'+j.id.replace('__10__','')] \
                                for j in rgenes) <=0, "knockout" + r.id + '_UB'
                    lp_prob += v[r.id] - (-M) * pulp.lpSum(1-z['u_G_'+j.id.replace('__10__','')] \
                                for j in rgenes) >=0, "knockout" + r.id + '_LB'
                # more complicated GPRs
                else:
                    # print r.id
                    # print GPR
                    proteins = [protein.replace("( ","").replace(" )","")\
                            .split(' and ') for protein in GPR.split(' or ')]
                    # print proteins
                    
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
                        gene = "u_G_" + gene.replace('__10__','')
                        lp_prob += v[r.id] - (1-z[gene])*M <= 0, \
                                    label + "_UB_" + r.id
                        lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, \
                                    label + "_LB_" + r.id
                        
                    allCombination = list(itertools.product(*nonCommonGenesList))
                    # print allCombination
                    # print allCombination
                    for i,genesC in enumerate(allCombination):
                        lp_prob += v[r.id] - M * pulp.lpSum(1-z['u_G_'+j.replace('__10__','')] \
                                    for j in genesC) <=0,\
                            "knockout" + r.id + '_UB_' + str(i)
                        lp_prob += v[r.id] - (-M) * pulp.lpSum(1-z['u_G_'+j.replace('__10__','')] \
                                    for j in genesC) >=0,\
                            "knockout" + r.id + '_LB_' + str(i)

        ############# constraints ################################
        print "add reaction indicator"    
        addReactionIndicator(lp_prob)

        #### ME model constraint
        S = self.get_S(me) # growth rate is 0.3
        # print S
        print "add GSM constraint"    
        # for i in metabolites:
        for i in S.keys():
            label = "mass_balance_%s"%i
            dot_S_v = pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()])
            condition = dot_S_v == 0
            lp_prob += condition, label

        ###### cut in the genome
        print "add cutting constraints"        
        lp_prob += pulp.lpSum(y[j] for j in genes) == 1, "end"
        lp_prob += pulp.lpSum(x[j] for j in genes) == 1, "start"

        # cut genes between start and end
        # for i,gene in enumerate(genes):
        #     lp_prob += pulp.lpSum(x[j] for j in \
        #               genes[0:i+1]) - pulp.lpSum(y[j] \
        #                   for j in genes[0:i+1]) - z[gene] == 0,\
        #                        'indicator' + str(gene)
        
        # A = pulp.LpAffineExpression()

        # for i,gene in enumerate(genes):
        #     A.addterm(x[gene],1) #pulp.lpSum(x[j] for j in genes[0:i+1])
        #     A.addterm(y[gene],-1)
        #     lp_prob += A - z[gene] == 0,'indicator' + str(gene)

        lp_prob += x[genes[0]] - y[genes[0]] == z[genes[0]], 'indicator' + str(genes[0])
        for i,gene in enumerate(genes):
            if i == 0: continue
            lp_prob += z[genes[i-1]] + x[gene] - y[gene] == z[gene],'indicator' + str(gene)


        ##### TUs 
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

        ##### some overlapped region cannot be selected as the start of deletion
        print "add no start and essential genes"
        for gene in no_start:
            lp_prob += x[gene] == 0, 'no_start_'+gene

        # # knock out transcription of cutted genes
        # for gene in genes:
        #     label = "knockout" + str(gene)
        #     # pdb.set_trace()
        #     if gene in v.keys():
        #         lp_prob += v[gene] - (1-z[gene])*M <= 0, label
        
        ##### add essential genes that cannot be deleted
        for eg in essential_genes:
            if eg in genes:
                lp_prob += z[eg] == 0
        
        ##### add regulation genes that cannot be deleted
        if reg_f != None:
            for eg in reg_genes:
                # remove the part joint with essential genes
                if (eg in genes) and (eg not in essential_genes):
                    lp_prob += z[eg] == 0

        ##### reaction bounds
        for r in me.reactions:
            # (lb,up) = me.bounds[r_id]
            v[r.id].lowBound = r.lower_bound
            v[r.id].upBound = r.upper_bound

        # r_biomass = 'bio00006'
        # v[r_biomass].lowBound = mu
        # lp_prob += v["bio00127_norm"] + v["bio00127b_norm"] >= mu 
        v['BIOMASS_Ec_iJO1366_WT_53p95M'].lowBound = mu

        v['BIOMASS_Ec_iJO1366_core_53p95M'].lowBound = 0
        v['BIOMASS_Ec_iJO1366_core_53p95M'].upBound = 0

        # lp file is somtime too larget to write
        # lp_prob.writeLP(lpfilename)

        # orignial implementation in the paper was calling cplex from C++ directly
        # call eternal compled cpp excutable to solve is a better option
        # it is implemented in codebase/mingenome_ecoli.cpp

        # current test version of using python to call the optimization
        # options = [epgap, epagap, epint, epopt, eprhs]
        GUROBI_CMD_OPTIONS = [('Threads', 8), ('TimeLimit', 1800), ('FeasibilityTol',1E-9),
                          ('OptimalityTol',1E-9),('IntFeasTol',1E-9),
                          ('MIPGapAbs', 0), ('MIPGap', 0), ('CliqueCuts', 2)]
        pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=0,
                                options=GUROBI_CMD_OPTIONS)
        # pulp_solver = pulp.solvers.CPLEX(path=None, keepFiles=0, mip=1,\
        #     msg=1, options=['mip tolerances mipgap 0', \
        #     'mip tolerances absmipgap 0', 'mip tolerances integrality 0',\
        #     'simplex tolerances optimality 1E-9',\
        #     'simplex tolerances feasibility 1E-9',], timelimit=1200)
        x_list = []
        y_list = []
        status = []
        def iterate_solve(lp_prob,iter_count):
            lp_prob.solve(pulp_solver)
            # print "----------- " + str(iter_count) + " ------------"
            status.append(pulp.LpStatus[lp_prob.status])
            # print("Status:", pulp.LpStatus[lp_prob.status])
            for v in lp_prob.variables():
                if "x_u_G_" in v.name and v.varValue == 1:
                    xname = v.name.replace("x_","")
                    # xname = xname.replace("J2_","J2-")
                    xname = xname.replace("PM0_","PM0-")
                    xname = xname.replace("PM_","PM-")
                    lp_prob += x[xname] == 1
                    if xname not in x_list: 
                        x_list.append(xname)
                if "y_u_G_" in v.name and v.varValue == 1:
                    yname = v.name.replace("y_","")
                    # yname = yname.replace("J2_","J2-")
                    yname = yname.replace("PM0_","PM0-")
                    yname = yname.replace("PM_","PM-")
                    lp_prob += y[yname] == 1
                    if yname not in y_list: 
                        y_list.append(yname)
            rhs = iter_count + 1
            lp_prob.constraints['start'].changeRHS(rhs)
            lp_prob.constraints['end'].changeRHS(rhs)
            return lp_prob
        
        for iter_count in xrange(1,12):
            lp_prob = iterate_solve(lp_prob,iter_count)
        pandas.DataFrame({'start': x_list, 'end':y_list,'status':status}).to_csv("./out/local_result_essential.csv")
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
    data_dir = "/data/SYN2973/2973.xml"
    model = modelCobra(data_dir)
    mip = build_MIP_by_Cobrapy(model, 0.98)

    #### analyze result
    # getAllgenesInResult()
