import vcf_utils.vcf_tools as vcf_tools
import math


def add_cnv_coverage(coverage_normalised, cnv_intervals, output_coverage):
    return True


def process_coverage(coverage_counts_file, coverage_bw_file, coverage_paired_normal, output_files):
    coverage_counts_handler = open(coverage_counts_file, 'r')
    chromosome_total, chromosome_list = find_chromosome_total(coverage_counts_handler)
    coverage_counts_handler.close()
    produce_normalised_coverage(coverage_bw_file, coverage_paired_normal, output_files)
    fill_coverage_circos_file(coverage_bw_file, output_files)


def get_position_coverage(counts_tumour, counts_normal, total_coverage_tumour, total_coverage_normal):
    rescaled_tumour_coverage = counts_tumour/total_coverage_tumour * 1000
    rescaled_normal_coverage = counts_normal/total_coverage_normal * 1000
    if (float(counts_normal)/100000) < 15:
        return 0
    if rescaled_normal_coverage == 0:
        return 0
    coverage_ratio = float(rescaled_tumour_coverage)/float(rescaled_normal_coverage)
    if coverage_ratio == 0:
        return 0
    coverage_log2 = math.log10(coverage_ratio)/math.log10(2)
    if coverage_log2 > 4:
        return 4
    return coverage_log2


def get_total_coverage(coverage_file):
    file = open(coverage_file, 'r')
    line = file.readline()
    total_counts = 0
    for line in file:
        line = line.split("\t")
        total_counts += (float(line[7]) * float(line[3]))
    file.close()
    return total_counts


def produce_normalised_coverage(coverage_bw_file, coverage_paired_normal, output_files):
    total_coverage_tumour = get_total_coverage(coverage_bw_file)
    total_coverage_normal = get_total_coverage(coverage_paired_normal)
    fileA = open(coverage_bw_file, 'r')
    fileA.readline()
    fileB = open(coverage_paired_normal, 'r')
    fileB.readline()
    objA = produce_bw_object(fileA.readline())
    objB = produce_bw_object(fileB.readline())
    chromosome = objA["chr"]
    start = 0
    fileB.readline()
    for lineA in fileA:
        objAA = produce_bw_object(lineA)
        lineB = fileB.readline()
        if lineB != "":
            objBB = produce_bw_object(lineB)
            normalised = get_position_coverage(objA["counts"], objB["counts"], total_coverage_tumour, total_coverage_normal)
            if chromosome == objAA["chr"]:
                if start + 100000 < objAA["end"]:
                    if normalised != 0:
                        output_files["coverage_normalised"].write(chromosome.replace("chr", "hs") + " " + str(start) + " " +
                                                              str(start + 100000) + " " +
                                                              str(normalised) + "\n")
                    start += 100000
                    objA = objAA
                    objB = objBB
                else:
                    objA["end"] = objAA["end"]
                    objA["counts"] += objAA["counts"]
                    objB["end"] = objBB["end"]
                    objB["counts"] += objBB["counts"]
            else:
                if normalised != 0:
                    output_files["coverage_normalised"].write(chromosome.replace("chr", "hs") + " " + str(start) + " " +
                                                          str(start + 100000) + " " + str(normalised) + "\n")
                chromosome = objAA["chr"]
                start = 0
                objA = objAA
                objB = objBB
    output_files["coverage_normalised"].close()


def produce_bw_object(line):
    line = line.split("\t")
    object = {}
    object["chr"] = line[0]
    object["start"] = line[1]
    object["end"] = line[2]
    object["size"] = line[3]
    object["mean"] = float(line[7])
    object["counts"] = float(object["mean"]) * float(object["size"])
    return object


def fill_coverage_circos_file(coverage_bw_file, output_files):
    coverage_bw_handler = open(coverage_bw_file, 'r')
    coverage_circos = output_files["coverage_distribution"]
    coverage_bw_handler.readline()
    basepairs = 0
    coverage_total = 0
    for line in coverage_bw_handler:
        line = line.split("\t")
        if float(line[7]) > 200:
            line[7] = 200
        if basepairs >= 100000:
            coverage_mean = coverage_total/100000
            coverage_circos.write(line[0].replace("chr", "hs") + " " + str(int(float(line[1]) - 50000)) +
                                  " " + str(int(float(line[1]) + 50000)) + " " + str(coverage_mean) + "\n")
            coverage_total = 0
            basepairs = 0
        else:
            basepairs = basepairs + float(line[2]) - float(line[1])
            coverage_total += (float(line[7]) * (float(line[2]) - float(line[1])))
    coverage_circos.close()
    coverage_bw_handler.close()


def produce_accumulated_coverage_histogram(coverage_file, chromosome_total, chromosome_list, coverage_histogram):
    chromosome_accum = initialize_chromosome_accumulated()
    coverage_file.readline()
    chromosome_accum = find_accumulated_frequencies(coverage_file, chromosome_accum, chromosome_total)
    chromosome_accum = process_and_clean_chromosome_accum(chromosome_accum, chromosome_list, coverage_histogram)
    return chromosome_accum


def process_and_clean_chromosome_accum(chromosome_accum, chromosome_list, coverage_histogram):
    t = 0
    while t < 24:
        i = 0
        while i < len(chromosome_accum[t]):
            chr_string = "chr" + str(t+1)
            chr_string = chr_string.replace("23", "X")
            chr_string = chr_string.replace("24", "Y")
            coverage_histogram.write(
                "hs" + str(chromosome_list[t]) + " " +
                str(vcf_tools.chromosome_length[chr_string]/300 * i) + " " +
                str(vcf_tools.chromosome_length[chr_string]/300 * (i + 1)) + " " +
                str(chromosome_accum[t+1][i]) + "\n"
            )
            i += 1
        t += 1
    return chromosome_accum


def find_accumulated_frequencies(coverage_file, chromosome_accum, chromosome_total):
    t = 0
    while t < 300:
        line = coverage_file.readline().split("\t")
        for i in range(0, len(chromosome_accum)):
            chromosome_accum[i].append(chromosome_accum[i][t] + (float(line[i])/chromosome_total[i] * 100))
        t += 1
    return chromosome_accum


def find_chromosome_total(coverage_file):
    chromosome_list = range(1, 23)
    chromosome_list.append("X")
    chromosome_list.append("Y")
    coverage_file.readline()
    line = coverage_file.readline()
    chromosome_total = map(float, line.split("\t"))
    for line in coverage_file:
        line = line.split("\t")
        i = 0
        while i < len(chromosome_total):
            chromosome_total[i] += float(line[i])
            i += 1
    coverage_file.close()
    return [chromosome_total, chromosome_list]


def initialize_chromosome_accumulated():
    t = 0
    chromosome_accum = []
    while t < 26:
        chromosome_accum.append([0])
        t += 1
    return chromosome_accum
