import os


def process_somatic_vcf(input_vcf, output_files, args):
    last = get_first_unfiltered_element(input_vcf)
    process_single_record(last, output_files)
    for record in input_vcf:
        if len(record.FILTER) == 0:
            process_single_record(record, output_files)
            process_consecutive_records(last, record, output_files, args.genome_version)
            last = record
    output_files["kataegis_out"].close()
    output_files["pre_indel_relative_density"].close()
    output_files["pre_snv_relative_density"].close()
    produce_relative_density(output_files["snv_relative_density"],
                             open(args.output_dir + "/pre_snv_relative_density.txt"))
    produce_relative_density(output_files["indel_relative_density"],
                             open(args.output_dir + "/pre_indel_relative_density.txt"))
    output_files["snv_positions"].close()
    snv_positions = open(args.output_dir + "/" + "snv_positions.txt", 'r')
    process_snv_context(snv_positions, output_files["snv_context"], args.genome_version)


def produce_relative_density(output_file, input_file):
    end = 2000000
    line = input_file.readline().split("\t")
    accumulated = 1
    chromosome = line[0]
    for line in input_file:
        line = line.split("\t")
        if line[0] == chromosome:
            if int(line[1]) > end and is_canonic_chromosome(chromosome):
                output_file.write(line[0].replace("chr", "hs") + " " + str(end - 1999999) + " " + str(end) + " " +
                                  str(accumulated) + "\n")
                end += 2000000
                accumulated = 1
            else:
                accumulated += 1
        else:
            if is_canonic_chromosome(chromosome):
                output_file.write(line[0].replace("chr", "hs") + " " + str(end - 2000000) + " " + str(end) + " " +
                                  str(accumulated) + "\n")
            accumulated = 1
            end = 2000000
            chromosome = line[0]
    output_file.close()


def get_first_unfiltered_element(input_vcf):
    last = next(input_vcf)
    while len(last.FILTER) > 0:
        last = next(input_vcf)
    return last


def process_structural_variation_vcf(input_vcf, output_files, args):
    for record in input_vcf:
        if len(record.FILTER) == 0:
            process_sv_record(record, output_files, args)
    output_files["traslocation_out"].close()
    output_files["inversion_out"].close()
    output_files["cnv_out"].close()
    output_files["duplication_out"].close()
    output_files["loh_out"].close()
    intersect_cnv_with_coverage(args, output_files)
    output_files["cnv_gain"].close()
    output_files["cnv_loss"].close()
    output_files["loh_region"].close()



def intersect_cnv_with_coverage(args, output_files):
    cnv_gain_file = open(args.output_dir + '/cnv_data.txt', 'r')
    for line in cnv_gain_file:
        line = line.split("\t")
        coverage_file = open(args.output_dir + '/coverage.txt', 'r')
        for cov in coverage_file:
            cov2 = cov.split(" ")
            if cov2[0].replace("hs", "chr") == line[0]:
                if cov2[1] >= line[1] and cov2[1] <= line[2]:
                    if int(line[3]) > 1:
                        output_files["cnv_gain"].write(cov)
                    else:
                        output_files["cnv_loss"].write(cov)
        coverage_file.close()
    loh_file = open(args.output_dir + '/loh_data.txt', 'r')
    for line in loh_file:
        line = line.split("\t")
        coverage_file = open(args.output_dir + '/coverage.txt', 'r')
        for cov in coverage_file:
            cov2 = cov.split(" ")
            if cov2[0].replace("hs", "chr") == line[0]:
                if cov2[1] >= line[1] and cov2[1] <= line[2]:
                   output_files["loh_region"].write(cov)
        coverage_file.close()


def process_sv_record(record, record_out, args):
    if "<INV" in str(record.ALT[0]):
        process_inversion(record, record_out)
    else:
        if "<CNV" in str(record.ALT[0]):
            process_copy_number(record, record_out)
        else:
            if "<DUP" in str(record.ALT[0]):
                process_duplication(record, record_out)
            else:
                if "<DEL" in str(record.ALT[0]):
                    process_deletion(record, record_out)
                else:
                    if "<INS" in str(record.ALT[0]):
                        return True
                    else:
                        if record.ALT == [None]:
                            return True
                        else:
                            process_structural_variant_manually(record, record_out)


def process_deletion(record, record_out):
    if abs(record.POS - record.INFO["END"]) >= 100000:
        record_out["deletion_out"].write("segdup" + str(record_out["segdup_counter"]) + " " + record.CHROM.replace("chr", "hs") + " " +
                                         str(record.POS) + " " + str(record.POS) + "\n")
        record_out["deletion_out"].write("segdup" + str(record_out["segdup_counter"]) + " " + record.CHROM.replace("chr", "hs") + " " +
                                         str(record.INFO["END"]) + " " + str(record.INFO["END"]) + "\n")
        record_out["segdup_counter"] += 1


def process_inversion(record, record_out):
    if abs(record.POS - record.INFO["END"]) >= 100000:
        record_out["inversion_out"].write("segdup" + str(record_out["segdup_counter"]) + " " +
                                          record.CHROM.replace("chr", "hs") + " " +
                                          str(record.POS) + " " + str(record.POS) + "\n")
        record_out["inversion_out"].write("segdup" + str(record_out["segdup_counter"]) + " " +
                                          record.CHROM.replace("chr", "hs") + " " +
                                          str(record.INFO["END"]) + " " + str(record.INFO["END"]) + "\n")
        record_out["segdup_counter"] += 1


def process_copy_number(record, record_out):
    if record.INFO["SVTYPE"] == "LOH":
        record_out["loh_out"].write(record.CHROM + "\t" + str(record.POS) + "\t" + str(record.INFO["END"]) + "\t" +
                                    str(record.samples[0]["CN"]) + "\n")
    else:
        record_out["cnv_out"].write(record.CHROM + "\t" + str(record.POS) + "\t" + str(record.INFO["END"]) + "\t" +
                                    str(record.samples[0]["CN"]) + "\n")
    record_out["segdup_counter"] += 1


def process_duplication(record, record_out):
    if abs(record.POS - record.INFO["END"]) >= 100000:
        record_out["duplication_out"].write("segdup" + str(record_out["segdup_counter"]) + " " +
                                            " ".join([record.CHROM.replace("chr", "hs"), str(record.POS),
                                            str(record.POS) + "\n"]) + "segdup" + str(record_out["segdup_counter"]) +
                                            " " + " ".join([record.CHROM.replace("chr", "hs"), str(record.INFO["END"]),
                                            str(record.INFO["END"]) + "\n"]))
        record_out["segdup_counter"] += 1


def happens_first(text, a, b):
    i = 0
    while i < len(text):
        if text[i] == a:
            return a
        if text[i] == b:
            return b
        i += 1
    return a


def get_allele_from_structural_variation(record):
    separator = happens_first(str(record.ALT[0]), "[", "]")
    return (str(record.ALT[0])).split(separator)[1]


def process_structural_variant_manually(record, record_out):
    return process_big_indel(record, record_out) if is_a_big_indel(record) else process_traslocation(record, record_out)


def is_a_big_indel(record):
    return all(map(is_nucleotide, record.REF)) and all(map(is_nucleotide, str(record.ALT[0])))


def is_canonic_chromosome(chromosome_string):
    return chromosome_string in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                                 "chr20", "chr21", "chr22", "chrX", "chrY"]


def process_traslocation(record, record_out):
    alternative = get_allele_from_structural_variation(record).split(":")
    if is_canonic_chromosome(record.CHROM) and is_canonic_chromosome(alternative[0]):
        alternative[0] = alternative[0].replace("chr", "hs")
        record_out["traslocation_out"].write(" ".join(["segdup" + str(record_out["segdup_counter"]),
                                                       record.CHROM.replace("chr", "hs"), str(record.POS),
                                                       str(record.POS) + "\n"]))
        record_out["traslocation_out"].write("segdup" + str(record_out["segdup_counter"]) + " " + alternative[0] + " " +
                                             str(alternative[1]) + " " + str(alternative[1]) + "\n")
        record_out["segdup_counter"] += 1


def process_big_indel(record, record_out):
    return True


def convert_to_absolute_genome_coordinates(chromosome, position, chromosome_accumulated):
    return chromosome_accumulated[chromosome] + position


def process_single_record(record, record_out):
    if is_snv(record):
        process_snv(record, record_out)
    else:
        if is_indel(record):
            process_indel(record, record_out)
        else:
            process_other(record, record_out)


def process_snv(record, record_out):
    record_out["snv_positions"].write(record.CHROM + "\t" + str(record.POS) + "\t" + record.REF + "\t" +
                                      str(record.ALT[0]) + "\n")
    record_out["snv_allele_frequencies"].write(str(find_snv_frequency(record)) + "\n")
    record_out["pre_snv_relative_density"].write(record.CHROM + "\t" + str(record.POS) + "\n")


def next_line_and_increment_counts(file_handler, counter):
    return [file_handler.readline().rstrip(), counter + 1]


def process_snv_context(snv_positions, output_file, genome_version):
    chromosome_file = open(os.path.dirname(os.path.realpath(__file__)) + "/../genome_files/" + genome_version +
                           "/chr1", 'r')
    chr_id = "chr1"
    chromosome_line, chromosome_line_count = next_line_and_increment_counts(chromosome_file, -1)
    chromosome_line_length = len(chromosome_line)
    counted = {}
    for line in snv_positions:
        line = line.rstrip().split("\t")
        pos_in_fasta_line = int(line[1]) % chromosome_line_length
        if is_canonic_chromosome(line[0]):
            if chromosomes_are_equal(line, chr_id):
                if pos_in_fasta_line == 1:
                    while chromosome_line_count < int(line[1])/chromosome_line_length:
                        prev_line = chromosome_line
                        chromosome_line, chromosome_line_count = next_line_and_increment_counts(chromosome_file,
                                                                                                chromosome_line_count)
                    counted = add_to_counted(line[2], line[3], prev_line[chromosome_line_length-1] +
                                             chromosome_line[0] + chromosome_line[1], counted)
                if pos_in_fasta_line == 0:
                    while chromosome_line_count < int(line[1])/chromosome_line_length:
                        prev_line = chromosome_line
                        chromosome_line, chromosome_line_count = next_line_and_increment_counts(chromosome_file,
                                                                                                chromosome_line_count)
                    counted = add_to_counted(line[2], line[3], prev_line[chromosome_line_length-2] +
                                           prev_line[chromosome_line_length - 1] + chromosome_line[0], counted)
                else:
                    while chromosome_line_count < int(line[1])/chromosome_line_length:
                        chromosome_line, chromosome_line_count = next_line_and_increment_counts(chromosome_file,
                                                                                                chromosome_line_count)
                    counted = add_to_counted(line[2], line[3], chromosome_line[pos_in_fasta_line - 2] +
                                             chromosome_line[pos_in_fasta_line-1] + chromosome_line[pos_in_fasta_line],
                                             counted)
            if chromosome_a_is_bigger(line, chr_id):
                chromosome_file.close()
                chromosome_file = open(os.path.dirname(os.path.realpath(__file__)) + "/../genome_files/" +
                                       genome_version + "/" + next_chromosome(chr_id), 'r')
                chr_id = next_chromosome(chr_id)
                chromosome_line, chromosome_line_count = next_line_and_increment_counts(chromosome_file, -1)
    mystring  = ""
    for change in ["A", "G", "T"]:
        for contc in ["ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT"]:
            mystring = mystring + str(counted[contc]["C"][change]) + "\t"
    for change in ["A", "C", "G"]:
        for contt in ["ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT"]:
            mystring = mystring + str(counted[contt]["T"][change]) + "\t"
    output_file.write(mystring + "\n")
    output_file.close()
    chromosome_file.close()


def add_to_counted(ref, alt, context, counted):
    if context[1] == "A" or context[1] == "G":
        context = complement_dna(reverse_seq(context))
        ref = complement_nucleotide(ref)
        alt = complement_nucleotide(alt)
    try:
        counted[context][ref][alt] += 1
    except:
        counted[context] = {"C": {"G": 0, "T": 0, "A": 0}, "T": {"A": 0, "G": 0, "C": 0}}
        counted[context][ref][alt] += 1
    return counted


def reverse_seq(seq):
    toreturn = ""
    i  = len(seq) -1
    while i >= 0:
        toreturn += seq[i]
        i -= 1
    return toreturn


def complement_nucleotide(nucl):
    complementary = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complementary[nucl]


def complement_dna(string):
    toreturn  = ""
    i = 0
    while i < len(string):
        toreturn += complement_nucleotide(string[i])
        i += 1
    return toreturn


def chromosome_a_is_bigger(line_a, line_b):
    return convert_chromosome_to_numeric(line_a[0]) > convert_chromosome_to_numeric(line_b)


def chromosomes_are_equal(line_a, line_b):
    return line_a[0] == line_b


def next_chromosome(chromosome):
    chr_id = convert_chromosome_to_numeric(chromosome)
    return convert_chromosome_from_numeric(chr_id + 1)


def convert_chromosome_to_numeric(chromosome):
    chromosomes = {"chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4, "chr5": 5, "chr6": 6, "chr7": 7, "chr8": 8, "chr9": 9,
                   "chr10": 10, "chr11": 11, "chr12": 12, "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16,
                   "chr17": 17, "chr18": 18, "chr19": 19, "chr20": 20, "chr21": 21, "chr22": 22, "chrX": 23,
                   "chrY": 24, "chrM": 25}
    return chromosomes[chromosome]


def convert_chromosome_from_numeric(chromosome):
    chromosomes = {1: "chr1", 2: "chr2", 3: "chr3", 4: "chr4", 5: "chr5", 6: "chr6", 7: "chr7", 8: "chr8", 9: "chr9",
                   10: "chr10", 11: "chr11", 12: "chr12", 13: "chr13", 14: "chr14", 15: "chr15", 16: "chr16",
                   17: "chr17", 18: "chr18", 19: "chr19", 20: "chr20", 21: "chr21", 22: "chr22", 23: "chrX",
                   24: "chrY", 25: "chrM"}
    return chromosomes[chromosome]


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def find_snv_frequency(record):
    return float(get_counts_for_alt_snv(record, 0))/(float(get_counts_for_alt_snv(record, 0)) +
                                                     float(get_counts_for_ref_snv(record, 0)))


def get_counts_for_ref_snv(record, sample_number):
    return record.samples[sample_number][str(record.REF) + "U"][0]


def get_counts_for_alt_snv(record, sample_number):
    return record.samples[sample_number][str(record.ALT[0]) + "U"][0]


def process_indel(record, record_out):
    record_out["indel_allele_frequencies"].write(str(find_indel_frequency(record)) + "\n")
    record_out["pre_indel_relative_density"].write(record.CHROM + "\t" + str(record.POS) + "\n")


def find_indel_frequency(record):
    return float(get_counts_for_alt_indel(record, 0)) / \
           (float(get_counts_for_alt_indel(record, 0)) + float(get_counts_for_ref_indel(record, 0)))


def get_counts_for_alt_indel(record, sample_number):
    return record.samples[sample_number]["TIR"][0]


def get_counts_for_ref_indel(record, sample_number):
    return record.samples[sample_number]["TAR"][0]


def process_other(record, record_out):
    record_out["other_out"].write(str(record) + "\n")


def process_consecutive_records(last, record, output_files, genome_version):
    if is_snv(last) and is_snv(record):
        chromosome_length, chromosome_accumulated = get_chromosome_length(genome_version)
        if last.CHROM in chromosome_length.keys():
            if find_distance(last, record) > 0:
                output_files["kataegis_out"].write(last.CHROM + "\t" +
                                                   str(convert_to_absolute_genome_coordinates(last.CHROM, last.POS,
                                                                                              chromosome_accumulated)) +
                                                   "\t" + str(find_distance(last, record)) + "\t" + str(last.REF) +
                                                   "\t" + str(last.ALT[0]) + "\t" + str(last.POS) + "\n")


def is_nucleotide(allele):
    return allele in ['A', 'T', 'C', 'G']


def are_all_nucleotides(record):
    return all(map(is_nucleotide, str(record.ALT[0]))) and all(map(is_nucleotide, str(record.REF)))


def lengths_dont_match(record):
    if record.ALT[0] is None:
        return 0 != len(record.REF)
    else:
        return len(record.ALT[0]) != len(record.REF)


def is_indel(record):
    return lengths_dont_match(record) and are_all_nucleotides(record)


def is_snv(record):
    return is_nucleotide(str(record.ALT[0])) and is_nucleotide(record.REF)


def find_distance(record_a, record_b):
    if record_a.CHROM == record_b.CHROM:
        return record_b.POS - record_a.POS
    return 0


def get_change(record_a):
    return str(record_a.REF[0]) + ">" + str(record_a.ALT[0])


def get_chromosome_length(genome_version):
    genome_size_file = open(os.path.dirname(os.path.realpath(__file__)) + "/../genome_files/" + genome_version +
                            "/chromosome_sizes", 'r')
    chromosome_size = {}
    chromosome_accumulated = {}
    total = 0
    for line in genome_size_file:
        line = line.split("\t")
        chromosome_size[line[0]] = line[1]
        chromosome_accumulated[line[0]] = total
        total += int(line[1])
    return [chromosome_size, chromosome_accumulated]
