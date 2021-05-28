#!/usr/bin/python
# -*-coding : utf-8-*-
# Copyright(c) 2018 - Fang Shuangsang <fangshuangsang@ict.ac.cn>
import math
import os
import sys
import optparse
import time
from multiprocessing import Pool
from sklearn.externals import joblib
import numpy as np
from collections import OrderedDict

sep = '\t'
alphabet = ['ttt', 'ttc', 'tta', 'ttg', 'tct', 'tcc', 'tca', 'tcg', 'tat', 'tac', 'taa', 'tag', 'tgt', 'tgc', 'tga',
            'tgg', 'ctt', 'ctc',
            'cta', 'ctg', 'cct', 'ccc', 'cca', 'ccg', 'cat', 'cac', 'caa', 'cag', 'cgt', 'cgc', 'cga', 'cgg', 'att',
            'atc', 'ata', 'atg', 'act', 'acc', 'aca', 'acg', 'aat', 'aac', 'aaa', 'aag', 'agt', 'agc', 'aga', 'agg',
            'gtt', 'gtc', 'gta', 'gtg', 'gct', 'gcc', 'gca', 'gcg', 'gat', 'gac', 'gaa', 'gag', 'ggt', 'ggc', 'gga',
            'ggg']


def cur_file_dir():
    path = sys.path[0]
    if os.path.isdir(path):
        return path
    elif os.path.isfile(path):
        return os.path.dirname(path)


curPath = cur_file_dir()


def cal_potential_score(trans_score, label):
    index_coding = '1'
    index_noncoding = '0'
    index = None
    # trans_score = trans_score * 1.0 / 30
    trans_score = math.log(trans_score, 10) if trans_score > 0 else 0.01
    if label == index_coding:
        index = 'coding'
        if trans_score < 0:
            trans_score = 0.01
    elif label == index_noncoding:
        index = 'noncoding'
        if trans_score > 0:
            trans_score = -1.0 / trans_score
        else:
            trans_score = -1.0 / 0.01
    potential_score = (sigmoid(trans_score) - 0.5) * 2
    # print trans_score, potential_score
    return index, potential_score


def sigmoid(strength):
    return 1 / (1 + np.exp(-strength))


def read_matrix(matrix_file):
    matrix_dict = {}
    with open(matrix_file, 'r') as fin:
        for line in fin:
            data = line.split()
            matrix_dict[data[0]] = float(data[1])
    return matrix_dict


def get_fasta_from_gtf(in_file, out_dir, two_bit):
    tmp_bed = '%s/in_file.bed' % out_dir
    out_fa = '%s/in_file.fa' % out_dir
    cmd1 = 'perl %s/gtf2Bed.pl %s > %s' % (curPath, in_file, tmp_bed)
    os.system(cmd1)
    cmd2 = '%s/twoBitToFa -bed=%s %s %s' % (curPath, tmp_bed, two_bit, out_fa)
    os.system(cmd2)
    return out_fa


def check_sequence(ori_id, ori_seq, log):
    new_seq = ori_seq.replace('u', 't').replace('\r', '')
    import re
    if re.search(r'[^atcg]', new_seq):
        log_str = ori_id + ' ' + 'contain unknow nucleotide, please checkout your sequence again. ' \
                                 'The transcipt will be passed\n'
        log.write(log_str)
        return False
    return True


def split_check_fasta_file(in_file, out_dir, parallel, log_file):
    two_line_fasta_list = []
    cur_id, cur_seq = None, ''
    log = open(log_file, 'w')
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.rstrip()
            if not line:
                break
            if line[0] == '>':
                if cur_id:
                    if check_sequence(cur_id, cur_seq, log):
                        two_line_fasta_list.append([cur_id, cur_seq])
                cur_id = line[1:]
                cur_seq = ''
            else:
                cur_seq = cur_seq + line.lower()
    if check_sequence(cur_id, cur_seq, log):
        two_line_fasta_list.append([cur_id, cur_seq])
    log.close()
    total_num = len(two_line_fasta_list)
    one_num = int(math.ceil(total_num * 1.0 / parallel))
    # print one_num
    out_file_prefix = '%s/CNCI2_file' % out_dir

    total_file_list = []
    for i in range(parallel):
        out_file = "%s_%s" % (out_file_prefix, str(i + 1))
        out = open(out_file, 'w')
        # print i * one_num, (i + 1) * one_num
        for data in two_line_fasta_list[i * one_num:(i + 1) * one_num]:
            out.write("%s\n" % sep.join(data))
        out.close()
        total_file_list.append(out_file)
    return total_file_list


def get_nt_list(seq):
    nt_list = []
    for i in range(len(seq) / 3):
        nt_list.append(seq[i * 3:i * 3 + 3])
    return nt_list


def get_nt_freq(nt_list):
    nt_freq = OrderedDict()
    for nt in alphabet:
        nt_freq[nt] = 0
    for nt in nt_list:
        nt_freq[nt] += 1
    total_nt_num = len(nt_list)

    nt_freq_list = []
    for k, v in nt_freq.items():
        nt_freq_list.append(str(v * 1.0 * len(alphabet) / total_nt_num))
    return nt_freq_list


def cal_score(seq, score_matrix, score_flag=False):
    coden_num = len(seq) / 3
    score = [-1000 for _ in range(coden_num)]
    start_list = [0]
    for i in range(1, coden_num):
        tmp_seq = seq[(i - 1) * 3: i * 3 + 3]
        tmp_score = score_matrix[tmp_seq]
        if score[i - 1] < 0:
            max_start = (i - 1)
            start_list.append(max_start)
            score[i] = tmp_score
        else:
            score[i] = score[i - 1] + tmp_score
            max_start = start_list[i - 1]
            start_list.append(max_start)
    # print score
    max_start, max_end = 0, 0
    max_score = -1000
    for i in range(len(score)):
        start = start_list[i] * 3
        end = i * 3 + 3
        cur_score = score[i]
        if cur_score > max_score:
            max_score = cur_score
            max_start, max_end = start + 1, end
    # print max_score, max_start, max_end
    if score_flag:
        return max_score, max_start, max_end, score
    return max_score, max_start, max_end


def get_feature(seq_score_list):
    nt_score = seq_score_list[0][0]
    std_score = np.std([float(score[0]) for score in seq_score_list])
    std_len = np.std([int(score[2]) - int(score[1]) +
                      1 for score in seq_score_list])
    nt_freq_list = [str(x) for x in seq_score_list[0][3]]
    return nt_score, std_score, std_len, nt_freq_list


def print_result(predict_file, detail_file, out_file):
    head = ['No.', 'Transcript ID', 'index', 'score', 'start', 'end', 'length']
    predict = open(predict_file, 'r')
    predict_list = [line.rstrip() for line in predict]
    predict.close()

    detail = open(detail_file, 'r')
    fout = open(out_file, 'w')
    fout.write("\t".join(head) + "\n")
    for i, line in enumerate(detail):
        data = line.rstrip().split(sep)
        trans_id, trans_start, trans_end, trans_score, trans_len = data
        trans_start, trans_end, trans_score, trans_len = int(trans_start), int(trans_end), float(trans_score), int(
            trans_len)
        # print trans_id, 'coding', str(trans_score), str(trans_start), str(trans_end), str(trans_len)
        # pre = trans_score
        # trans_score = trans_score / (trans_len / 3 - 1)
        # print pre, trans_score, trans_len / 3
        index, potential_score = cal_potential_score(
            trans_score, predict_list[i])
        trans_num = str(i+1)
        fout.write("\t".join(
            [trans_num, trans_id, index, str(potential_score), str(trans_start), str(trans_end), str(trans_len)]) + "\n")
    fout.close()


def main(in_file, matrix_dict, out_dir, score_flag=False):
    temp_feature_file = os.path.join(out_dir, 'CNCI2_score')
    temp_detail_file = os.path.join(out_dir, 'CNCI2_detail')
    feature_out = open(temp_feature_file, 'w')
    detail_out = open(temp_detail_file, 'w')
    if score_flag:
        mlcds_score_out = os.path.join(out_dir, 'mlcds_score_out')
        if not os.path.exists(mlcds_score_out):
            os.mkdir(mlcds_score_out)
    with open(in_file, 'r') as fin:
        cnt = 0
        for line in fin:
            if not line:
                break
            cnt += 1
            trans_id, trans_seq = line.rstrip().split(sep)
            seq_len = len(trans_seq)
            seq_score_list = []
            for pos in range(6):
                if pos <= 2:
                    new_seq = trans_seq[pos: (seq_len - pos) / 3 * 3 + pos]
                else:
                    pos = pos - 3
                    new_seq = trans_seq[::-
                                        1][pos: (seq_len - pos) / 3 * 3 + pos]
                # start, end are real position based on 1
                res = cal_score(new_seq, matrix_dict, score_flag=score_flag)
                # print res
                if score_flag:
                    score, start, end, score_list = res
                else:
                    score, start, end = res
                    score_list = None
                nt_list = get_nt_list(new_seq)
                nt_freq_list = get_nt_freq(nt_list)
                seq_score_list.append(
                    (score, start + pos, end + pos, nt_freq_list, score_list))
            seq_score_list = sorted(
                seq_score_list, key=lambda s: s[0], reverse=True)

            mlcds = seq_score_list[0]
            [max_score, max_start, max_end] = mlcds[:3]
            detail_str = "%s\n" % sep.join([str(trans_id), str(max_start), str(max_end), str(
                max_score), str(len(trans_seq))])
            detail_out.write(detail_str)

            nt_score, std_score, std_len, nt_freq_list = get_feature(
                seq_score_list)
            feature_list = [
                str(x) for x in [nt_score, std_score, std_len, "\t".join(nt_freq_list)]]
            feature_str = "%s\n" % sep.join(feature_list)
            feature_out.write(feature_str)

            if score_flag:
                tmp_id = trans_id.split()[0]
                mlcds_out_file = os.path.join(
                    mlcds_score_out, '%s.score.csv' % cnt)
                mlcds_data_tmp = [tmp[4] for tmp in seq_score_list]
                min_len = min([len(x) for x in mlcds_data_tmp])
                mlcds_data_tmp = [x[:min_len] for x in mlcds_data_tmp]
                # mlcds_data = np.array(mlcds_data_tmp)
                with open(mlcds_out_file, 'w') as fout:
                    fout.write(
                        ",".join(['mlcds', 'np1', 'np2', 'np3', 'np4', 'np5']) + "\n")
                    for i in range(1, min_len):
                        # print len(mlcds_data_tmp)
                        tmp = ','.join([str(x[i]) for x in mlcds_data_tmp])
                        fout.write("%s\n" % tmp)
                fout.close()
                mlcds_out_pic = os.path.join(
                    mlcds_score_out, '%s.score.png' % cnt)
                r_plot_code = os.path.join(curPath, 'plot.R')
                os.system("Rscript %s %s %s" %
                          (r_plot_code, mlcds_out_file, mlcds_out_pic))
                # print mlcds_data[:10]
                # # mlcds_data = mlcds_data.T
                # print mlcds_data[:10]
                # exit()
                # with open(mlcds_out_file, 'w') as fout:
                #     for tmp in seq_score_list:
                #         score_list = tmp[4]
                #         fout.write("%s\n" % ",".join([str(x) for x in score_list[1:]]))
                # fout.close()
                # head = ",".join(['mlcds', 'np1', 'np2', 'np3', 'np4', 'np5'])
                # print 'save'
                # np.savetxt(mlcds_out_file, mlcds_data, delimiter=',',  header=head)

        feature_out.close()
        detail_out.close()
        return temp_feature_file, temp_detail_file


def merge_file(score_file, detail_file, index_total_res):
    total_score = open(score_file, 'w')
    total_detail = open(detail_file, 'w')
    for (tmp_score_file, tmp_detail_file) in index_total_res:
        with open(tmp_score_file, 'r') as tmp_score:
            for line in tmp_score:
                total_score.write(line)
        tmp_score.close()

        with open(tmp_detail_file, 'r') as tmp_detail:
            for line in tmp_detail:
                total_detail.write(line)
            tmp_detail.close()
    total_score.close()
    total_detail.close()


def run():
    in_tmp_file = options.file
    out_dir = options.outfile
    parallel = 1
    class_model = options.model
    gtf_type = options.gtf
    genome_two_bit = options.genome_two_bit
    score_flag = False
    remove_tmp_dir = True

    start_time = time.time()
    tmp_out_dir = os.path.join(out_dir, 'Temp_Dir')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists(tmp_out_dir):
        os.mkdir(tmp_out_dir)
    log_file = os.path.join(out_dir, 'CNCI2.log')
    print "CNCI2 start index ..."
    if gtf_type:
        if not genome_two_bit:
            exit("when gtf file input, genome two bit is needed.")
        in_file = get_fasta_from_gtf(in_tmp_file, tmp_out_dir, genome_two_bit)
    else:
        in_file = in_tmp_file

    matrix_file, model = None, None
    if class_model == 've':
        matrix_file = os.path.join(curPath, "CNCI2_Parameters/animal.matrix")
        model = os.path.join(curPath, "CNCI2_Parameters/animal.model")
    elif class_model == 'pl':
        matrix_file = os.path.join(curPath, "CNCI2_Parameters/plant.matrix")
        model = os.path.join(curPath, "CNCI2_Parameters/plant.model")
    else:
        exit("wrong model was selected. ve(vertebrate) or pl(plant)")

    matrix_dict = read_matrix(matrix_file)
    total_file_list = split_check_fasta_file(
        in_file, tmp_out_dir, parallel, log_file)
    score_file, detail_file = None, None
    if parallel == 1:
        in_file = total_file_list[0]
        if not os.path.getsize(in_file):
            print "No sequences to predict!"
            exit()
        score_file, detail_file = main(
            in_file, matrix_dict, tmp_out_dir, score_flag)
    elif parallel > 1:
        index_pool = Pool(parallel)
        index_result = []
        for i, in_file in enumerate(total_file_list):
            index_result.append(index_pool.apply_async(
                main, args=(in_file, matrix_dict, out_dir, score_flag)))
        index_pool.close()
        index_pool.join()
        index_total_res = [res.get() for res in index_result]
        score_file = os.path.join(tmp_out_dir, 'CNCI2_score')
        detail_file = os.path.join(tmp_out_dir, 'CNCI2_detail')
        merge_file(score_file, detail_file, index_total_res)

    run_model = joblib.load(model)
    tmp_score = np.loadtxt(score_file, dtype=np.float)
    if len(tmp_score.shape) == 1:
        tmp_score = np.array([tmp_score])
    predict_out = run_model.predict(tmp_score)
    predict_file = "%s_pred" % in_file
    np.savetxt(predict_file, predict_out, fmt='%d')
    final_result_file = os.path.join(out_dir, "%s.index" % out_dir)
    print_result(predict_file, detail_file, final_result_file)
    if remove_tmp_dir:
        import shutil
        shutil.rmtree(tmp_out_dir, True)
    end_time = time.time()
    print "CNCI2 index completed. Elapsed: %ss " % (end_time - start_time)


if __name__ == "__main__":
    print("CNCIv2 2018/11/16\n")
    parse = optparse.OptionParser()
    parse.add_option('-f', '--file', dest='file', action='store', metavar='input files',
                     help='enter your transcript (sequence or gtf)')
    parse.add_option('-o', '--out', dest='outfile', action='store', metavar='output files',
                     help='assign your output file')
    # parse.add_option('-p', '--parallel', dest='parallel', action='store', metavar='prallel numbers',
    #                  help='please enter your specified speed ratio')
    parse.add_option('-m', '--model', dest='model', action='store', metavar='model types', default='ve',
                     help='please enter your specified classification model')
    parse.add_option('-g', '--gtf', dest='gtf', action='store_true', metavar='gtf file name',
                     help='please enter your gtf files')
    parse.add_option('-b', '--genome_two_bit', dest='genome_two_bit', action='store', metavar='genome two bit file',
                     help='if your input file is gtf type please enter RefGenome two bit file')
    (options, args) = parse.parse_args()
    run()
