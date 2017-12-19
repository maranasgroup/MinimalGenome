from Bio import SeqIO
import pdb
import pandas
import json

class utility(object):
    """docstring for utility"""
    def __init__(self, sequence_dir):
        super(utility, self).__init__()
        self.gb_file = sequence_dir
    
    def getAllGenesfromGB(self):
        for record in SeqIO.parse(self.gb_file, "genbank"):
            print "locus_tag", "old_locus_tag", "gene", "start", "end", "strand"
            for f in record.features:
                if f.type == "gene":
                    locus_tag = f.qualifiers["locus_tag"][0]
                    if "old_locus_tag" not in f.qualifiers:
                        old_locus_tag = "NA"
                    else:
                        old_locus_tag = f.qualifiers["old_locus_tag"][0]

                    if "gene" not in f.qualifiers:
                        gene = "NA"
                    else:
                        gene = f.qualifiers["gene"][0]

                    start = f.location.start.position
                    end = f.location.end.position
                    strand= f.strand

                    print locus_tag, old_locus_tag, gene, start+1, end, strand

    
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

    def move_start_site_to_back(self,data_dir):
        all = data_dir + "/genes_and_promoters.csv"
        df = pandas.read_csv(all)

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
        pandas.DataFrame(start_list).to_csv(data_dir + "/genes_and_promoters_new_start.csv")

    def exchange_promoter_key_value(self,data_dir):
        all = data_dir + "/genes_and_promoters.xlsx"
        df = pandas.read_excel(all)
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
        with open(data_dir + '/gene_promoter_dict.json', 'w') as fp:
            json.dump(gene_promoter_dict, fp)

    def get_promoter_end_loc(self):
        all = "/Users/linuswang/Documents/MinimalGenome/src/data/Bacillus/genes_and_promoters.xlsx"
        df = pandas.read_excel(all)
        dict_all = df.to_dict(orient='index')

        for key in dict_all:
            if dict_all[key]["class"] == "promoter":
                # search forward
                # pdb.set_trace()
                if dict_all[key]["strand"] == 1:
                    key_tmp = key
                    while dict_all[key_tmp]["class"] != "gene":
                        key_tmp = key_tmp + 1
                    dict_all[key]['end'] = dict_all[key_tmp]['start'] - 1
                
                # search backward
                if dict_all[key]["strand"] == -1:
                    key_tmp = key
                    while dict_all[key_tmp]["class"] != "gene":
                        key_tmp = key_tmp - 1
                    dict_all[key]['end'] = dict_all[key]['start']
                    dict_all[key]['start'] = dict_all[key_tmp]['end'] + 1                  

        # print dict_all[0]
        df_new = pandas.DataFrame.from_dict(dict_all,orient='index')
        df_new.to_csv("/Users/linuswang/Documents/MinimalGenome/src/data/Bacillus/genes_and_promoters.csv")

    def map_genes_between_models(self,data_dir):
        all = data_dir + "/map_genes.xlsx"
        df_1 = pandas.read_excel(all,sheetname="Sheet1")
        df_2 = pandas.read_excel(all,sheetname="sub_wiki")
        genes_alias = df_1['GENE ALIASES'].tolist()
        all_alias = df_2['other_names'].tolist()
        locus_tags = df_2['locus_tag'].tolist()
        # print genes_alias
        indexes = []

        for gene in genes_alias:
            index = 0
            success = 0
            for genes in all_alias:
                index = index + 1
                if gene in genes:
                    success = 1
                    break
            if success == 0: 
                index = 99999999
            indexes.append(index)
            print gene, index
    
    def draw(self,file_name,sheetname):
        df = pandas.read_excel(file_name,sheetname=sheetname)
        starts = df["start_loc"].tolist()
        ends = df["end_loc"].tolist()
        print starts
        print ends
        for i,start in enumerate(starts):
            if i == 25: break
            print '<trackmarker start="' + str(start) + '" end="'+ str(ends[i]) + \
        '" markerstyle="fill:rgba(170,85,0,0.9)" length="4" startlength="-4"></trackmarker>'

if __name__ == '__main__':
    # find_genes_in_TU()
    # get_overlap()
    # remove_duplicate_promoter()
    # get_promoter_removed_df()
    # move_start_site_to_back()
    # exchange_promoter_key_value()
    
    # #### solve
    # data_dir = "../data/SYN2973/2973.xml"
    # model = modelCobra(data_dir)
    # mip = build_MIP_by_Cobrapy(model, 0.98)

    #### analyze result
    # getAllgenesInResult()

    # draw deletions
    file_name = "../../doc/result_deletions_bacillus.xlsx"
    sheetname="model_normalized_daniel"#"minGenome_normalized"# "experimental_PG10" # "exprimental_PG38"
    u = utility('')
    u.draw(file_name,sheetname)