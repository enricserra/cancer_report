import os
import subprocess


def prepare_genome_files(genome_version, genome_ref):
    if exists_genome_files(genome_version):
        return True
    else:
        print "Needed genome files for " + genome_version + " not found, creating them \n"
        subprocess.call("mkdir   " + get_genome_files_dir() + "/" + genome_version, shell=True)
        parse_genome(genome_version, genome_ref)


def exists_genome_files(genome_version):
    if os.path.isdir(get_genome_files_dir() + genome_version):
        return True
    else:
        return False


def parse_genome(genome_version, genome_ref):
    genome = open(genome_ref, 'r')
    chromosome_sizes = open(get_genome_files_dir() + genome_version +
                            "/chromosome_sizes", 'w')
    chromosome_size = 0
    previous_chromosome = ""
    for genome_line in genome:
        if is_chromosome(genome_line):
            if is_canonic_chromosome(previous_chromosome):
                chromosome_sizes.write(previous_chromosome + "\t" + str(chromosome_size) + "\n")
            chromosome_name = get_chromosome_name(genome_line)
            print "Creating files for " + chromosome_name
            previous_chromosome = chromosome_name
            chromosome_size = 0
            try:
                chromosome_out.close()
                chromosome_out = open(get_genome_files_dir() + genome_version + "/" + chromosome_name, 'w')
            except NameError:
                chromosome_out = open(get_genome_files_dir() + genome_version + "/" + chromosome_name, 'w')
        else:
            chromosome_out.write(genome_line)
            chromosome_size += (len(genome_line) -1)


def get_genome_files_dir():
    return os.path.dirname(os.path.realpath(__file__)) + "/../genome_files/"


def is_chromosome(line):
    return line[0] == ">"


def get_chromosome_name(genome_line):
    return genome_line.strip(">").rstrip()


def is_canonic_chromosome(chromosome_string):
    return chromosome_string in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                                 "chr20", "chr21", "chr22", "chrX", "chrY"]



