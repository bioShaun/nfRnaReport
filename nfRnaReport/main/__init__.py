# coding:UTF-8
import os
import jinja2
from copy import deepcopy
import configparser
'''
init html and pdf configure
'''
# template env config
html_jinja_env = jinja2.Environment(
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader(os.path.join(
        os.path.dirname(__file__), 'html_templates'))
)
pdf_jinja_env = jinja2.Environment(
    block_start_string='\BLOCK{',
    block_end_string='}',
    variable_start_string='\VAR{',
    variable_end_string='}',
    comment_start_string='\#{',
    comment_end_string='}',
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader(os.path.join(
        os.path.dirname(__file__), 'pdf_templates'))
)

##################
# html version
##################

CFG_DIR = os.path.dirname(os.path.realpath(__file__))

# add each analysis part
configFilePath = os.path.join(os.path.dirname(__file__), 'report_conf.conf')
command = configparser.ConfigParser()
command.read(configFilePath)
# all path
mRNA_data_path = command.get('mRNA-report-result', 'report_result_path')
mRNA_result_path = command.get('mRNA-report-result', 'report_result_path')
enrichment_path = command.get('mRNA-path', 'enrichment_path')
fastqc_path = command.get('mRNA-path', 'fastqc_path')
mapping_path = command.get('mRNA-path', 'mapping_path')
quantification_path = command.get('mRNA-path', 'quantification_path')
rseqc_path = command.get('mRNA-path', 'rseqc_path')
nr_assembly_path = command.get('mRNA-path', 'nr_assembly_path')
novel_tr_path = command.get('mRNA-path', 'novel_tr_path')
transfactor_path = command.get('mRNA-path', 'transfactor_path')
rmats_path = command.get('mRNA-path', 'rmats_path')
snp_path = command.get('mRNA-path', 'snp_path')
ppi_path = command.get('mRNA-path', 'ppi_path')


# pdf settings
address = command.get('mRNA-pdf-info', 'address')
phone = command.get('mRNA-pdf-info', 'phone')
table_rows = command.get('mRNA-pdf-info', 'table_rows')
max_cell_len = command.get('mRNA-pdf-info', 'max_cell_len')
project_name = command.get('mRNA-pdf-info', 'project_name')
logo_path = os.path.join(CFG_DIR, command.get(
    'mRNA-pdf-static', 'logo_path'))
cover_path = os.path.join(CFG_DIR, command.get(
    'mRNA-pdf-static', 'cover_path'))
pipeline_path = os.path.join(CFG_DIR, command.get(
    'mRNA-pdf-static', 'pipeline_path'))
mRNAworkflow_path = os.path.join(CFG_DIR, command.get(
    'mRNA-pdf-static', 'mRNAworkflow_path'))
rmaths_path = os.path.join(CFG_DIR, command.get(
    'mRNA-pdf-static', 'rmats'
))
rmaths_inclusion_path = os.path.join(CFG_DIR, command.get(
    'mRNA-pdf-static', 'rmaths_inclusion_path'
))


mRNA_data_dict = dict(enrichment=enrichment_path,
                      fastqc=fastqc_path,
                      mapping=mapping_path,
                      quantification=quantification_path,
                      rseqc=rseqc_path,
                      nr_assembly=nr_assembly_path,
                      novel_tr=novel_tr_path,
                      transfactor=transfactor_path,
                      rmats=rmats_path,
                      snp=snp_path,
                      ppi=ppi_path)

for key, value in mRNA_data_dict.items():
    mRNA_data_dict[key] = os.path.join(mRNA_data_path, value)

mRNA_result_dict = deepcopy(mRNA_data_dict)
for key, value in mRNA_result_dict.items():
    mRNA_result_dict[key] = os.path.join(mRNA_result_path, value)

#################
# pdf version
#################

# enrichment part
enrichment_analysis_path = dict(
    go_barplot_path='go/*/*go.enrichment.barplot.png',
    kegg_barplot_path='kegg/*/*kegg.enrichment.barplot.png',
    pathview_path='kegg/*/*pathway/src/*.png',
    dag_bp_path='go/*/DAG/*BP.GO.DAG.png',
    dag_cc_path='go/*/DAG/*CC.GO.DAG.png',
    dag_mf_path='go/*/DAG/*MF.GO.DAG.png',
    go_table_path='go/*/*go.enrichment.txt',
    kegg_table_path='kegg/*/*kegg.enrichment.txt')

for key, value in enrichment_analysis_path.items():
    enrichment_analysis_path[key] = os.path.join(
        mRNA_data_path, enrichment_path, value)

# fastqc part
fastqc_analysis_path = dict(
    gc_barplot_path='reads_gc/gc_distribution.line.report.png',
    reads_quality_path='reads_quality/reads_quality.bar.report.png',
    qc_table_path='data.summary.csv',
    reads_filter_path='reads_filter/reads_filter.pie.report.png')

for key, value in fastqc_analysis_path.items():
    fastqc_analysis_path[key] = os.path.join(
        mRNA_data_path, fastqc_path, value)

# mapping part
mapping_analysis_path = dict(
    mapping_table_path='overall_mapping_stats.txt',
    mapping_plot_path='mapping_portion.png')

for key, value in mapping_analysis_path.items():
    mapping_analysis_path[key] = os.path.join(
        mRNA_data_path, mapping_path, value)

# quantification part
quantification_analysis_path = dict(
    correlation_heatmap_path='expression_summary/Sample.correlation.heatmap.png',
    gene_expression_path='expression_summary/Gene_expression.png',
    pca_plot_path='expression_summary/PCA_plot.png',
    gene_table_path='expression_summary/Gene.tpm.txt',
)

for key, value in quantification_analysis_path.items():
    quantification_analysis_path[key] = os.path.join(
        mRNA_data_path, quantification_path,
        value)

# diff part
diff_analysis_path = dict(
    volcano_plot_path='differential_analysis/*/*Volcano_plot.png',
    diff_heatmap_path='expression_summary/Diff.genes.heatmap.png',
    diff_table_path='differential_analysis/*/*.edgeR.DE_results.txt')

for key, value in diff_analysis_path.items():
    diff_analysis_path[key] = os.path.join(
        mRNA_data_path, quantification_path, value)
# rseqc part
rseqc_analysis_path = dict(
    genebody_coverage_plot_path='genebody_coverage/genebody_coverage.report.png',
    seq_saturation_plot_path='sequencing_saturation/sequence_saturation.report.png')

for key, value in rseqc_analysis_path.items():
    rseqc_analysis_path[key] = os.path.join(mRNA_data_path, rseqc_path, value)

# nr_assembly part
nr_assembly_analysis_path = dict(
    iso_len_dis_path='Isoform_length_distribution.png',
    gene_len_dis_path='Gene_length_distribution.png',
    nr_assembly_stats='Trinity.stat.txt',
)

for key, value in nr_assembly_analysis_path.items():
    nr_assembly_analysis_path[key] = os.path.join(
        mRNA_data_path, nr_assembly_path, value)

# novel transcript
novel_tr_analysis_path = dict(
    novel_tr_table_path='novel.gene.annotation.csv'
)
for key, value in novel_tr_analysis_path.items():
    novel_tr_analysis_path[key] = os.path.join(
        mRNA_data_path, novel_tr_path, value)

# transfactor
transfactor_analysis_path = dict(
    transfactor_number_table_path='Diff_TF_number.stats.csv',
    transfactor_detail_table_path='TF_in_diff_compare.stats.csv'
)
for key, value in transfactor_analysis_path.items():
    transfactor_analysis_path[key] = os.path.join(
        mRNA_data_path, transfactor_path, value)


# rmats
rmats_analysis_path = dict(
    rmats_bar_plot='*/*as_summary.png',
)
for key, value in rmats_analysis_path.items():
    rmats_analysis_path[key] = os.path.join(
        mRNA_data_path, rmats_path, value)

# snp
snp_analysis_path = dict(
    snp_table_path='all_sample.vcf.table.txt'
)
for key, value in snp_analysis_path.items():
    snp_analysis_path[key] = os.path.join(
        mRNA_data_path, snp_path, value)

# snp
ppi_analysis_path = dict(
    ppi_table_path='*.ppi.csv'
)
for key, value in ppi_analysis_path.items():
    ppi_analysis_path[key] = os.path.join(
        mRNA_data_path, ppi_path, value)


pdf_analysis_path = {'enrichment': enrichment_analysis_path,
                     'fastqc': fastqc_analysis_path,
                     'mapping': mapping_analysis_path,
                     'quantification': quantification_analysis_path,
                     'diff': diff_analysis_path,
                     'rseqc': rseqc_analysis_path,
                     'nr_assembly': nr_assembly_analysis_path,
                     'novel_tr': novel_tr_analysis_path,
                     'transfactor': transfactor_analysis_path,
                     'rmats': rmats_analysis_path,
                     'snp': snp_analysis_path,
                     'ppi': ppi_analysis_path}


# other info
pdf_settings = {'address': address,
                'phone': phone,
                'table_rows': table_rows,
                'max_cell_len': max_cell_len,
                'logo_path': logo_path,
                'cover_path': cover_path,
                'project_name': project_name,
                'pipeline_path': pipeline_path,
                'mRNAworkflow_path': mRNAworkflow_path,
                'rmaths_path': rmaths_path,
                'rmaths_inclusion_path': rmaths_inclusion_path}
