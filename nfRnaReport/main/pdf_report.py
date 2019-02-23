# coding:UTF-8
'''
this is py2report's pdf_report moudle which a python script generate pdf mRNA report
this version only run on 34
log:
create by chencheng on 2017-06-13
add all href and plot_size on 2017-06-14
add all function doc on 2017-06-16
'''
import os
import sys
import subprocess
from . import mRNA_data_dict, pdf_analysis_path
from . import pdf_jinja_env, pdf_settings
import glob

reload(sys)
sys.setdefaultencoding('utf-8')

main_path = os.path.dirname(os.path.realpath(__file__))
ref_file_path = os.path.join(main_path, 'pdf_templates', 'ref.bib')


def cut_overlong_table(row_list, max_len=int(pdf_settings['max_cell_len'])):
    '''
    param:
    row_list: each analysis part's table rwo list apply in three_line_list function
    max_cell_len: threshold value to cut table cell
    function: cut table cell's length when over 20 chars
    '''
    for i in range(len(row_list)):
        if len(row_list[i]) > max_len:
            row_list[i] = row_list[i][:max_len] + '...'
    return row_list


def three_line_list(input_path,
                    split='\t', limited_rows=True,
                    colunms=None,
                    max_len=int(pdf_settings['max_cell_len'])):
    '''
    param:
    input_path: each analysis part's table path
    colunms: columns to include in table, number start from 0
    function: transform each analysis part's table into table_list
    '''
    if input_path:
        with open(input_path, 'r+') as f:
            data = f.readlines()
            thead = data[0]
            table_cols = len(thead.strip().split(split))
            tbody = data[1:]
            if len(tbody) > int(pdf_settings['table_rows']) and limited_rows:
                tbody = tbody[:int(pdf_settings['table_rows'])]

            if colunms is None:
                cols_num = table_cols
            else:
                cols_num = len(colunms)

            table_list = []
            table_begin = '\\begin{tabular}{%s}' % ('c' * cols_num)
            table_list.append(table_begin)

            head_list = cut_overlong_table(
                thead.strip('\n').split(split))
            if colunms is not None:
                head_list = [head_list[each] for each in colunms]
            head_str = '&'.join(head_list).replace('_', '\_').replace(
                '%', '\%').replace('#', '\#') + r'\\'
            table_list.append(head_str)
            for line in tbody:
                each_list = cut_overlong_table(
                    line.strip('\n').split(split),
                    max_len=max_len)
                if colunms is not None:
                    each_list = [each_list[each] for each in colunms]
                each_str = '&'.join(each_list).replace('_', '\_').replace(
                    '%', '\%').replace('#', '\#') + r'\\'
                table_list.append(each_str)
            return table_list
    else:
        return None


def check_file(file_dict, generate_report_path, part):
    '''
    param:
    file_dict:each analysis part's path dict
    generate_report_path:a path where to your analysis's report data
    function:check each path's exists and add generate_report_path into them
    '''
    for key, value in file_dict.items():
        file_dict[key] = os.path.join(generate_report_path, value)
        file_list = glob.glob(file_dict[key])
        if not file_list or not os.path.exists(file_list[0]):
            if part:
                print '{file} is not find in: {file_path}'.format(
                    file=key, file_path=file_list[0])
                sys.exit(1)
            else:
                file_dict[key] = None
        else:
            file_dict[key] = file_list[0]


def run_tex(tex_path):
    '''
    param:
    tex_path:generate pdf report's tex path
    function:transform tex file into pdf file
    '''
    tex_dir = os.path.dirname(tex_path)
    tex_file = os.path.basename(tex_path)
    os.chdir(tex_dir)
    aux_file = tex_file.replace('tex', 'aux')
    rm_set = ['.aux', '.log', '.out', '.toc',
              '.tmp', '.bib', '.bbl', '.blg', '.tex']
    subprocess.call('cp {ref_file} {tex_dir}'.format(
        ref_file=ref_file_path, tex_dir=tex_dir), shell=True)
    subprocess.call('xelatex {tex_file} > summary'.format(
        tex_file=tex_file), shell=True)
    subprocess.call('bibtex {aux_file}'.format(aux_file=aux_file), shell=True)
    subprocess.call('xelatex {tex_file} > summary'.format(
        tex_file=tex_file), shell=True)
    subprocess.call('xelatex {tex_file} > summary'.format(
        tex_file=tex_file), shell=True)
    for each_file in os.listdir(tex_dir):
        if os.path.splitext(each_file)[1] in rm_set:
            subprocess.call('rm {file}'.format(
                file=os.path.join(tex_dir, each_file)), shell=True)

    print '---------------------'
    print 'pdf mRNA report done!'
    print '---------------------'


def enrichment_analysis_part(generate_report_path, part):
    enrichment_analysis_path = pdf_analysis_path['enrichment']
    check_file(enrichment_analysis_path, generate_report_path, part)
    kegg_list = three_line_list(
        enrichment_analysis_path['kegg_table_path'], colunms=range(7))
    go_list = three_line_list(
        enrichment_analysis_path['go_table_path'], colunms=range(7))

    if (enrichment_analysis_path['dag_bp_path'] and
        enrichment_analysis_path['dag_cc_path'] and
            enrichment_analysis_path['dag_mf_path']):
        dag_plots = True
    else:
        dag_plots = False

    enrichment_dict = dict(
        enrichment='enrichment', kegg_begin=kegg_list[0],
        kegg_head=kegg_list[1], kegg_body=kegg_list[2:],
        go_begin=go_list[0], go_head=go_list[1], go_body=go_list[2:],
        go_barplot_path=enrichment_analysis_path['go_barplot_path'],
        dag_plots=dag_plots,
        dag_bp_path=enrichment_analysis_path['dag_bp_path'],
        dag_cc_path=enrichment_analysis_path['dag_cc_path'],
        dag_mf_path=enrichment_analysis_path['dag_mf_path'],
        kegg_barplot_path=enrichment_analysis_path['kegg_barplot_path'],
        pathview_path=enrichment_analysis_path['pathview_path'],
        go_table_path=enrichment_analysis_path['go_table_path'],
        kegg_table_path=enrichment_analysis_path['kegg_table_path'])

    return enrichment_dict


def fastqc_analysis_part(generate_report_path, part):
    fastqc_analysis_path = pdf_analysis_path['fastqc']
    check_file(fastqc_analysis_path, generate_report_path, part)
    qc_list = three_line_list(fastqc_analysis_path['qc_table_path'], colunms=range(7),
                              limited_rows=False, split=',')

    fastqc_dict = dict(
        qc_begin=qc_list[0], qc_head=qc_list[1], qc_body=qc_list[2:],
        gc_barplot_path=fastqc_analysis_path['gc_barplot_path'],
        reads_quality_path=fastqc_analysis_path['reads_quality_path'],
        qc_table_path=fastqc_analysis_path['qc_table_path'],
        data_stat='data_stat',
        reads_filter_path=fastqc_analysis_path['reads_filter_path']
    )

    return fastqc_dict


def mapping_analysis_part(generate_report_path, part):
    mapping_analysis_path = pdf_analysis_path['mapping']
    check_file(mapping_analysis_path, generate_report_path, part)
    mapping_list = three_line_list(
        mapping_analysis_path['mapping_table_path'], colunms=range(5))

    mapping_dict = dict(
        mapping_begin=mapping_list[0],
        mapping_head=mapping_list[1],
        mapping_body=mapping_list[2:],
        mapping_plot_path=mapping_analysis_path['mapping_plot_path'],
        mapping_table_path=mapping_analysis_path['mapping_table_path'],
        mapping='mapping')

    return mapping_dict


def rseqc_analysis_part(generate_report_path, part):
    rseqc_analysis_path = pdf_analysis_path['rseqc']
    if os.path.exists(
            os.path.join(generate_report_path, mRNA_data_dict['rseqc'])):
        check_file(rseqc_analysis_path, generate_report_path, part)
        rseqc_dict = dict(
            data_control='data_control',
            genebody_coverage_plot_path=rseqc_analysis_path['genebody_coverage_plot_path'],
            seq_saturation_plot_path=rseqc_analysis_path['seq_saturation_plot_path'],
        )

        return rseqc_dict


def novel_tr_analysis_part(generate_report_path, part):
    novel_tr_analysis_path = pdf_analysis_path['novel_tr']
    if os.path.exists(
        os.path.join(generate_report_path, mRNA_data_dict['novel_tr'])
    ):
        check_file(novel_tr_analysis_path, generate_report_path, part)
        novel_tr_cols = range(6)
        novel_tr_list = three_line_list(
            novel_tr_analysis_path['novel_tr_table_path'],
            colunms=novel_tr_cols,
            split=',')
        novel_tr_dict = dict(
            novel_tr_begin=novel_tr_list[0],
            novel_tr_head=novel_tr_list[1],
            novel_tr_body=novel_tr_list[2:],
            novel_tr='novel_tr',
            novel_tr_table_path=novel_tr_analysis_path['novel_tr_table_path'])
        return novel_tr_dict


def transfactor_analysis_part(generate_report_path, part):
    transfactor_analysis_path = pdf_analysis_path['transfactor']
    if os.path.exists(
        os.path.join(generate_report_path, mRNA_data_dict['transfactor'])
    ):
        check_file(transfactor_analysis_path, generate_report_path, part)
        tr_num_list = three_line_list(
            transfactor_analysis_path['transfactor_number_table_path'],
            split=',')
        tr_detail_list = three_line_list(
            transfactor_analysis_path['transfactor_detail_table_path'],
            split=',',
            colunms=range(7)
        )
        transfactor_dict = dict(
            tr_num_begin=tr_num_list[0],
            tr_num_head=tr_num_list[1],
            tr_num_body=tr_num_list[2:],
            tr_detail_begin=tr_detail_list[0],
            tr_detail_head=tr_detail_list[1],
            tr_detail_body=tr_detail_list[2:16],
            transfactor='transfactor',
            transfactor_number_table_path=transfactor_analysis_path['transfactor_number_table_path'],
            transfactor_detail_table_path=transfactor_analysis_path['transfactor_detail_table_path'])
        return transfactor_dict


def rmats_analysis_part(generate_report_path, part):
    rmats_analysis_path = pdf_analysis_path['rmats']
    check_file(rmats_analysis_path, generate_report_path, part)
    rmats_dict = dict(
        rmats_bar_plot=rmats_analysis_path['rmats_bar_plot'],
        rmats='rmats')

    return rmats_dict


def snp_analysis_part(generate_report_path, part):
    snp_analysis_path = pdf_analysis_path['snp']
    if os.path.exists(
        os.path.join(generate_report_path, mRNA_data_dict['snp'])
    ):
        check_file(snp_analysis_path, generate_report_path, part)
        snp_list = three_line_list(
            snp_analysis_path['snp_table_path'],
            colunms=[0, 1, 2, 3, 6, 8, 9])
        snp_dict = dict(
            snp_begin=snp_list[0],
            snp_head=snp_list[1],
            snp_body=snp_list[2:],
            snp='snp',
            snp_table_path=snp_analysis_path['snp_table_path'])
        return snp_dict


def ppi_analysis_part(generate_report_path, part):
    ppi_analysis_path = pdf_analysis_path['ppi']
    if os.path.exists(
        os.path.join(generate_report_path, mRNA_data_dict['ppi'])
    ):
        check_file(ppi_analysis_path, generate_report_path, part)
        ppi_list = three_line_list(
            ppi_analysis_path['ppi_table_path'],
            colunms=[0, 1, 2, 3, 17],
            split=',')
        ppi_dict = dict(
            ppi_begin=ppi_list[0],
            ppi_head=ppi_list[1],
            ppi_body=ppi_list[2:],
            ppi='ppi',
            ppi_table_path=ppi_analysis_path['ppi_table_path'])
        return ppi_dict


def nr_assembly_analysis_part(generate_report_path, part):
    nr_assembly_path = pdf_analysis_path['nr_assembly']
    if os.path.exists(
            os.path.join(generate_report_path, mRNA_data_dict['nr_assembly'])):
        check_file(nr_assembly_path, generate_report_path, part)
        nr_assembly_s_list = three_line_list(
            nr_assembly_path['nr_assembly_stats'], colunms=3
        )
        nr_assembly_dict = dict(
            nr_assembly_stats_begin=nr_assembly_s_list[0],
            nr_assembly_stats_head=nr_assembly_s_list[1],
            nr_assembly_stats_body=nr_assembly_s_list[2:],
            nr_assembly_stats=nr_assembly_path['nr_assembly_stats'],
            iso_len_dis_path=nr_assembly_path['iso_len_dis_path'],
            gene_len_dis_path=nr_assembly_path['gene_len_dis_path'],
            nr_assembly='nr_assembly',
        )

        return nr_assembly_dict


def quant_analysis_part(generate_report_path, part):
    diff_analysis_path = pdf_analysis_path['diff']
    quantification_analysis_path = pdf_analysis_path['quantification']
    check_file(diff_analysis_path, generate_report_path, part)
    check_file(quantification_analysis_path, generate_report_path, part)
    diff_cols = [0, -3, -2, -1]
    diff_list = three_line_list(
        diff_analysis_path['diff_table_path'], colunms=diff_cols)
    gene_count_list = three_line_list(
        quantification_analysis_path['gene_table_path'], colunms=range(6))
    quant_dict = dict(
        gene_count_begin=gene_count_list[0],
        gene_count_head=gene_count_list[1],
        gene_count_body=gene_count_list[2:],
        diff_begin=diff_list[0],
        diff_head=diff_list[1],
        diff_body=diff_list[2:],
        correlation_heatmap_path=quantification_analysis_path['correlation_heatmap_path'],
        gene_expression_path=quantification_analysis_path['gene_expression_path'],
        pca_plot_path=quantification_analysis_path['pca_plot_path'],
        gene_table_path=quantification_analysis_path['gene_table_path'],
        volcano_plot_path=diff_analysis_path['volcano_plot_path'],
        diff_heatmap_path=diff_analysis_path['diff_heatmap_path'],
        diff_table_path=diff_analysis_path['diff_table_path'],
        quant='quant',
        diff='diff')

    return quant_dict


def check_analysis_part(generate_report_path, analysis_part,
                        part_dict, label, part, func):
    if os.path.exists(os.path.join(generate_report_path, analysis_part)):
        analysis_dict = func(generate_report_path, part)
    else:
        analysis_dict = part_dict
        print '---{analysis_module} analysis part missing---'.format(
            analysis_module=label)

    return analysis_dict
    # pdf_param_dict.update(analysis_dict)


def create_pdf_report(generate_report_path, part):
    '''
    param:a path where to your analysis's report data
    function:generate report tex file
    '''
    pdf_param_dict = {}
    pdf_project_name = generate_report_path.rstrip(
        '/').rsplit('/', 1)[1].replace('_', '\_')
    pdf_head_dict = dict(
        project_name=pdf_project_name,
        report_name=pdf_settings['project_name'],
        address=pdf_settings['address'],
        phone=pdf_settings['phone'],
        logo_path=pdf_settings['logo_path'],
        cover_path=pdf_settings['cover_path'],
        pipeline_path=pdf_settings['pipeline_path'],
        mRNAworkflow_path=pdf_settings['mRNAworkflow_path'],
        rmaths_path=pdf_settings['rmaths_path'],
        rmaths_inclusion_path=pdf_settings['rmaths_inclusion_path']
    )
    pdf_param_dict.update(pdf_head_dict)
    # enrichment:
    enrichment_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['enrichment'],
        {'enrichment': None},
        'enrichment', part,
        func=enrichment_analysis_part)
    pdf_param_dict.update(enrichment_dict)
    # fastqc:
    fastqc_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['fastqc'],
        {'data_stat': None},
        'fastqc', part,
        func=fastqc_analysis_part)
    pdf_param_dict.update(fastqc_dict)
    # mapping:
    mapping_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['mapping'],
        {'mapping': None},
        'mapping', part,
        func=mapping_analysis_part)
    pdf_param_dict.update(mapping_dict)
    # rseqc:
    rseqc_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['rseqc'],
        {'data_control': None},
        'rseqc', part,
        func=rseqc_analysis_part)
    pdf_param_dict.update(rseqc_dict)
    # nr_assembly
    nr_assembly_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['nr_assembly'],
        {'nr_assembly': None},
        'nr_assembly',
        part,
        func=nr_assembly_analysis_part
    )
    pdf_param_dict.update(nr_assembly_dict)
    # novel tr
    novel_tr_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['novel_tr'],
        {'novel_tr': None},
        'novel_tr',
        part,
        func=novel_tr_analysis_part
    )
    pdf_param_dict.update(novel_tr_dict)
    # transfactor
    transfactor_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['transfactor'],
        {'transfactor': None},
        'transfactor',
        part,
        func=transfactor_analysis_part
    )
    pdf_param_dict.update(transfactor_dict)
    # rmats
    rmats_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['rmats'],
        {'rmats': None},
        'rmats',
        part,
        func=rmats_analysis_part
    )
    pdf_param_dict.update(rmats_dict)
    # snp
    snp_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['snp'],
        {'snp': None},
        'snp',
        part,
        func=snp_analysis_part
    )
    pdf_param_dict.update(snp_dict)
    # snp
    ppi_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['ppi'],
        {'ppi': None},
        'ppi',
        part,
        func=ppi_analysis_part
    )
    pdf_param_dict.update(ppi_dict)
    # quant&diff
    quant_dict = check_analysis_part(
        generate_report_path,
        mRNA_data_dict['quantification'],
        {'quant': None, 'diff': None},
        'quant&diff', part,
        func=quant_analysis_part)
    pdf_param_dict.update(quant_dict)

    template = pdf_jinja_env.get_template('mRNA_base')
    pdf_template_path = generate_report_path

    if not os.path.exists(pdf_template_path):
        os.makedirs(pdf_template_path)
    with open(os.path.join(pdf_template_path,
                           'report.tex'), 'w+') as f:
        f.write(template.render(pdf_param_dict))
    print '-------------------------'
    print 'pdf report tex file done!'
    print '-------------------------'
    run_tex(tex_path=os.path.join(pdf_template_path, 'report.tex'))
