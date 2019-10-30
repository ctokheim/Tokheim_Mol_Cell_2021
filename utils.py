import pandas as pd
import numpy as np
import scipy.stats as stats
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

#############################
# Categorize UPS function
#############################

def categorize_ups(gene_list):
    # read in all the annotation
    e1 = pd.read_csv('data/literature_ups/E1.txt', sep='\t').Gene.unique()
    e2 = pd.read_csv('data/literature_ups/E2.txt', sep='\t').Gene.unique()
    e3 = pd.read_csv('data/literature_ups/E3.txt', sep='\t').Gene.unique()
    dub = pd.read_csv('data/literature_ups/DUB.txt', sep='\t').Gene.unique()
    target_recognition = pd.read_csv('data/literature_ups/Target_recognition.txt', sep='\t').Gene.unique()
    
    # mark a dataframe for each driver's categories
    ups_df = pd.DataFrame({'gene': list(gene_list), 
                           'E3_Ge_et_al': 0, 
                           'Target recognition': 0,
                           'E2': 0,
                           'E1': 0,
                           'DUB': 0})
    ups_df.loc[ups_df['gene'].isin(e3), 'E3_Ge_et_al'] = 1
    ups_df.loc[ups_df['gene'].isin(target_recognition), 'Target recognition'] = 1
    ups_df.loc[ups_df['gene'].isin(e2), 'E2'] = 1
    ups_df.loc[ups_df['gene'].isin(e1), 'E1'] = 1
    ups_df.loc[ups_df['gene'].isin(dub), 'DUB'] = 1
    ups_df['E3'] = ups_df[['E3_Ge_et_al', 'Target recognition']].max(axis=1)
    return ups_df

#############################
# UPS driver gene functions
#############################

def read_driver_gene_results(path='data/ups_driver_gene/2020plus/'):
    output_dict = {}
    mypattern = os.path.join(path, '*.txt')
    for f in glob.glob(mypattern):
        # read in data
        tmp_df = pd.read_csv(f, sep='\t')
        
        # filter data
        is_signif = tmp_df['qvalue']<=0.05
        is_high_effect = tmp_df['score']>0.5
        signif_genes = tmp_df.loc[is_signif & is_high_effect, 'gene'].unique()
        
        # save results in dictionary
        ctype = os.path.basename(f)[:-4]
        output_dict[ctype] = signif_genes
    return output_dict

def read_all_results_ups_dataframe(base_dir='data/ups_driver_gene/2020plus/',
                                   qval_thresh=0.05):
    # add cancer types
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir))]
    output_list = []
    for c in cancer_types:
        tmp = pd.read_table('{0}/{1}.txt'.format(base_dir, c))
        tmp['cancer type'] = c
        tmp = tmp[(tmp['qvalue']<=qval_thresh) & (tmp['score']>0.5)]
        output_list.append(tmp)
    # merge results
    merged_df = pd.concat(output_list)
    merged_df['gene type'] = merged_df['info'].str.extract('TYPE=(.+)$')
    return merged_df

def merge_dictionary(dict1, dict2):
    """Merge the result of two dictionaries"""
    all_keys = set(dict1.keys()) | set(dict2.keys())
    output_dict = {}
    for k in all_keys:
        output_dict[k] = set(dict1.get(k, [])) | set(dict2.get(k, []))
    return output_dict

def read_all_results_ub_sites(base_dir='data/degron/ub_site_enrichment',
                              qval_thresh=0.1,
                              col_name='qvalue'):
    """Reads all of the results"""
    # add cancer types
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir))]

    # append results for various cancer types
    output_list = []
    for c in cancer_types:
        tmp = pd.read_table('{0}/{1}.txt'.format(base_dir, c))
        tmp['cancer type'] = c
        tmp = tmp[tmp[col_name]<=qval_thresh]
        output_list.append(tmp)

    # merge results
    merged_df = pd.concat(output_list)

    return merged_df

def read_all_results_phosphodegron_sites(base_dir='data/degron/phosphodegron_enrichment',
                              qval_thresh=0.1,
                              col_name='qvalue'):
    """Reads all of the results"""
    # add cancer types
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir))]

    # append results for various cancer types
    output_list = []
    for c in cancer_types:
        tmp = pd.read_table('{0}/{1}.txt'.format(base_dir, c))
        tmp['cancer type'] = c
        tmp = tmp[tmp[col_name]<=qval_thresh]
        output_list.append(tmp)

    # merge results
    merged_df = pd.concat(output_list)

    return merged_df

def read_all_results_known_degron(base_dir='data/degron/known_degron_enrichment',
                                  qval_thresh=0.1):
    """Reads all of the results"""
    # add cancer types
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir))]

    # append results for various cancer types
    output_list = []
    for c in cancer_types:
        tmp = pd.read_table('{0}/{1}.txt'.format(base_dir, c))
        tmp['cancer type'] = c
        tmp = tmp[tmp['qvalue']<=qval_thresh]
        output_list.append(tmp)

    # merge results
    merged_df = pd.concat(output_list)

    return merged_df

def read_mutation_flags(path='data/degron/mutation_flags_degron/'):
    output_list = []
    mypattern = os.path.join(path, '*.txt')
    for f in glob.glob(mypattern):
        ctype = os.path.basename(f)[:-4]
        tmp = pd.read_table(f, index_col=0)
        tmp_sum = tmp.sum()
        for gene in tmp_sum.index:
            output_list.append([gene, ctype, tmp_sum.loc[gene]])
    mut_ct_df = pd.DataFrame(output_list, columns=['gene', 'cancer type', 'count'])
    return mut_ct_df

###############
# Venn diagram
###############

def overlap_set(set1, set2, n=20177):
    # overlap sets
    len_intersect = len(set(set1) & set(set2))
    len_set1 = len(set1)
    len_set2 = len(set2)
    
    # return contigency table
    contig_tbl = [
        [len_intersect, len_set1 - len_intersect],
        [len_set2 - len_intersect, n-len(set(set1) | set(set2))]
    ]
    return contig_tbl

def venn_diagram(set1, set2, name1, name2, title=''):
    len_intersect = len(set(set1) & set(set2))
    len_set1 = len(set1)
    len_set2 = len(set2)
    overlap = (len_set1 - len_intersect, len_set2 - len_intersect, len_intersect)

    # plot venn diagram
    with sns.plotting_context('paper', font_scale=2.5):
        venn2(subsets=overlap, set_labels=(name1, name2))
        venn2_circles(subsets=overlap, linestyle='solid', linewidth=.75)
        plt.title(title, size=26)
        
    return overlap

def venn_diagram3(set1, set2, set3, name1, name2, name3, title='', ax=None):
    # make sure to convert to set object
    set1 = set(set1)
    set2 = set(set2)
    set3 = set(set3)
    
    # do all possible intersections
    intersect_12 = set(set1) & set(set2)
    intersect_13 = set(set1) & set(set3)
    intersect_23 = set(set2) & set(set3)
    full_intersect = set(set1) & set(set2) & set(set3)
    
    # figure out number only specific to one set
    len_set1_specific = len(set1 - set2 - set3)
    len_set2_specific = len(set2 - set1 - set3)
    len_set3_specific = len(set3 - set1 - set2)
    
    # Figure out length of full intersect
    len_full_intersect = len(full_intersect)
    
    # figure out length of two-set specific overlaps
    len_set12 = len(intersect_12 - full_intersect)
    len_set13 = len(intersect_13 - full_intersect)
    len_set23 = len(intersect_23 - full_intersect)

    # create overlap object
    #overlap = (len_set1 - len_intersect, len_set2 - len_intersect, len_intersect)
    overlap = {'100': len_set1_specific, '010': len_set2_specific, '001': len_set3_specific,
               '110': len_set12, '101': len_set13, '011': len_set23,
               '111': len_full_intersect}

    # plot venn diagram
    with sns.plotting_context('notebook', font_scale=1.0):
        if ax:
            venn3(subsets=overlap, set_labels=(name1, name2, name3), ax=ax)
            venn3_circles(subsets=overlap, linestyle='solid', linewidth=.75, ax=ax)
        else:
            venn3(subsets=overlap, set_labels=(name1, name2, name3))
            venn3_circles(subsets=overlap, linestyle='solid', linewidth=.75)            
        plt.title(title, size=16)
        
    return overlap
        
#########################
# Read CGC
#########################

def process_cgc(path, return_dataframe=False):
    """Get the list of CGC genes with small somatic variants."""
    # read in data
    df = pd.read_table(path)

    # keep small somatic variants
    s = df['Mutation Types']
    is_small = s.str.contains('Mis|F|N|S').fillna(False)
    is_somatic = ~df['Tumour Types(Somatic)'].isnull()
    df = df[is_small & is_somatic].copy()

    # get gene names
    if not return_dataframe:
        cgc_genes = df['Gene Symbol'].tolist()
    else:
        cgc_genes = df

    return cgc_genes

#############################
# GPS assay
#############################
from sklearn.model_selection import train_test_split

def bulky(mystring):
    mysum = 0
    for x in mystring:
        if x in ['F', 'W', 'Y']:
            mysum += 1
    return mysum

def acidic(mystring):
    mysum = 0
    for x in mystring:
        if x in ['D', 'E']:
            mysum += 1
    return mysum

def compare_motif_cterm(peptide, motif):
    for i in range(len(motif)):
        if peptide[-(i+1)]!=motif[-(i+1)] and motif[-(i+1)]!='x':
            return 0
    return 1

def compare_motif_nterm(peptide, motif):
    for i in range(len(motif)):
        if peptide[i]!=motif[i] and motif[i]!='x':
            return 0
    return 1

def motif_count(mystring, motif_list, motif_type='cterm'):
    mysum = 0
    for motif in motif_list:
        if motif_type == 'cterm':
            mysum += compare_motif_cterm(mystring, motif)
        else:
            mysum += compare_motif_nterm(mystring, motif)
    return mysum

def rule_based(data_path='data/GPS/cterm_deepDegron/cterm_degron_predictions.txt', 
               top100_path='data/GPS/cterm_degron_koren_et_al_rule_based.txt',
               motif_type='cterm'):
    """Add simple rule-based approaches to predictions of GPS assay"""
    # read in the data
    df = pd.read_csv(data_path, sep='\t')
    
    # get top 100 motifs
    rule_degron_df = pd.read_csv(top100_path, sep='\t')
    top100_motifs = rule_degron_df['Motif'].iloc[:100].tolist()

    # rule-based predictors
    df['bulky'] = df['Peptide amino acid sequence'].apply(bulky)
    df['acidic'] = -df['Peptide amino acid sequence'].apply(acidic)
    df['simple motif'] = df['Peptide amino acid sequence'].apply(motif_count, args=(top100_motifs, motif_type))
    df['y'] = (df['Modal Bin']<2.5).astype(int)

    # recapitulate split of training/test
    train_x, test_x, train_y, test_y = train_test_split(df, df['y'], 
                                                        train_size=0.7, 
                                                        test_size=0.3, 
                                                        random_state=101, 
                                                        shuffle=True)
    
    return train_x, test_x, train_y, test_y

def process_df(mydf):
    mydf['alt'] = mydf['Substitution'].str[-1:]
    mydf['ref'] = mydf['Substitution'].str[:1]
    mydf['pos'] = mydf['Substitution'].str[1:-1].astype(int)
    baseline = mydf[mydf['ref']==mydf['alt']]['degron potential'].iloc[0]
    mydf['delta degron potential'] = mydf['degron potential'] - baseline
    return mydf

def create_sat_mutagenesis_summary(path):
    import glob ; import os
    prng = np.random.RandomState(101)
    output_list = []
    mypattern = os.path.join(path, '*.txt')
    for f in glob.glob(mypattern):
        bootstrap_list = []
        gene = os.path.basename(f)[:-4]
        tmp_df = pd.read_csv(f, sep='\t')
        tmp_df = process_df(tmp_df)
        r, pval = stats.pearsonr(tmp_df['delta degron potential'], tmp_df['PSI'])

        # bootstrap to get confidence intervals
        for i in range(1000):
            boot_df = tmp_df.sample(frac=1, replace=True, random_state=prng)
            r_boot, pval_boot = stats.pearsonr(boot_df['delta degron potential'], boot_df['PSI'])
            bootstrap_list.append(r_boot)
        myseries = pd.Series(bootstrap_list)
        lower = myseries.quantile(0.025)
        upper = myseries.quantile(0.975)

        # append to output
        output_list.append([gene, r, lower, upper, pval])
    
    # create summary output
    summary_df = pd.DataFrame(output_list, columns=['gene', 'pearson r', 'lower', 'upper', 'pval'])
    summary_df.sort_values('pearson r', inplace=True)
    
    return summary_df

#######################
# immune-related funcs
#######################
def read_immune_data(mut_ct_df, immune_path, count_thresh=5, qval_thresh=0.1):
    # read in immune association
    immune_df = pd.read_csv(immune_path, sep='\t')
    split = immune_df.analysis.str.split('_', expand=True).rename(columns={0: 'cancer type', 1: 'gene'})
    immune_df = pd.concat([immune_df, split], axis=1)
    immune_df = pd.merge(immune_df, mut_ct_df, on=['gene', 'cancer type'], how='left')

    # perform analysis of immune-related biomarkers
    lf_df = immune_df[(immune_df['variable']=='LeukocyteFraction') & (immune_df['count']>count_thresh)].copy()
    lf_df['qvalue'] = bh_fdr(lf_df['pvalue'])
    tgfb_df = immune_df[(immune_df['variable']=='TGFB') & (immune_df['count']>count_thresh)].copy()
    tgfb_df['qvalue'] = bh_fdr(tgfb_df['pvalue'])
    ifng_df = immune_df[(immune_df['variable']=='IFNG') & (immune_df['count']>count_thresh)].copy()
    ifng_df['qvalue'] = bh_fdr(ifng_df['pvalue'])
    macrophage_df = immune_df[(immune_df['variable']=='MacrophageRegulation') & (immune_df['count']>count_thresh)].copy()
    macrophage_df['qvalue'] = bh_fdr(macrophage_df['pvalue'])
    wound_df = immune_df[(immune_df['variable']=='WoundHealing') & (immune_df['count']>count_thresh)].copy()
    wound_df['qvalue'] = bh_fdr(wound_df['pvalue'])
    lymph_df = immune_df[(immune_df['variable']=='LymphocyteInfiltration') & (immune_df['count']>count_thresh)].copy()
    lymph_df['qvalue'] = bh_fdr(lymph_df['pvalue'])
    leuk_df = immune_df[(immune_df['variable']=='LeukocyteFraction') & (immune_df['count']>count_thresh)].copy()
    leuk_df['qvalue'] = bh_fdr(leuk_df['pvalue'])
    
    # filter to only significant
    signif_immune = pd.concat([
         lymph_df[lymph_df['qvalue']<qval_thresh],
         #leuk_df[leuk_df['qvalue']<qval_thresh],
         tgfb_df[tgfb_df['qvalue']<qval_thresh],
         ifng_df[ifng_df['qvalue']<qval_thresh],
         macrophage_df[macrophage_df['qvalue']<qval_thresh],
         wound_df[wound_df['qvalue']<qval_thresh]
    ])
    uniq_signif = signif_immune['analysis'].unique()
    variable_list = ['LymphocyteInfiltration', 
                     'TGFB', 'IFNG', 'MacrophageRegulation', 'WoundHealing']
    
    # create heatmap dataframe
    is_variable = immune_df['variable'].isin(variable_list)
    is_signif = immune_df['analysis'].isin(uniq_signif)
    plot_df = immune_df[is_variable & is_signif].pivot(index='variable', columns='analysis', values='tvalue')
    
    return signif_immune, plot_df


########################
# stat funcs
########################
def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.
    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.
    Parameters
    ----------
    pval : list or array
        list/array of p-values
    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]
    
    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(int(n))
    i = np.arange(1, int(n)+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]