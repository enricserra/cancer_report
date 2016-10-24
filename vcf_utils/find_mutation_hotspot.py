import subprocess


def process_kataegis(kataegis_input_file, output_files, r_scripts, output_dir):
    kataegis_chromosomes, distance_median = generate_kataegis_chromosome_and_median(kataegis_input_file)
    hotspots = process_kataegis_chromosomes(kataegis_chromosomes, distance_median)
    output_files["hotspots_json"].write("[")
    for chromosome in hotspots.keys():
        for pos in hotspots[chromosome]:
            output_files["hotspots_file"].write(kataegis_chromosomes[chromosome][pos].split("\t")[1] + "\n")
            output_files["hotspots_json"].write("\"" + str(chromosome) + ":" +
                                                str(convert_absolute_position_to_relative(chromosome, kataegis_chromosomes[chromosome][pos - 5].split("\t")[1])) + "," + 
                                                str(convert_absolute_position_to_relative(chromosome, kataegis_chromosomes[chromosome][pos - 1].split("\t")[1])) +
                                                 "\", ")
    output_files["hotspots_file"].close()
    output_files["hotspots_json"].write("]")
    output_files["hotspots_json"].close()
    call_r_plotting_scripts(r_scripts, output_dir)


def convert_absolute_position_to_relative(chromosome, absolute_position):
    chromosome_accumulated = [0, 248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912,
                           1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562,
                           2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417, 3088269832]
    chromosome = chromosome.replace("chr", "")
    if chromosome == "X":
        chromosome = 23
    if chromosome == "Y":
        chromosome = 24
    if chromosome == "M":
        chromosome = 25
    return int(absolute_position) - int(chromosome_accumulated[int(chromosome) - 1])


def call_r_plotting_scripts(r_scripts, output_dir):
    subprocess.call("Rscript " + r_scripts["plot_frequencies"] + " " + output_dir + "/snv_frequencies.txt " +
                    output_dir + "/indel_frequencies.txt " + " " + output_dir + "/", shell=True)
    subprocess.call("Rscript " + r_scripts["plot_kataegis"] + " " + output_dir + "/kata.txt " +
                    output_dir + "/kataegis_hotspots.txt " + output_dir, shell=True)


def generate_kataegis_chromosome_and_median(kataegis_input_file):
    median_list = []
    kataegis_chromosomes = {}
    for line in kataegis_input_file:
        chromosome = line.split("\t")[0]
        median_list.append(int(line.split("\t")[2]))
        try:
            kataegis_chromosomes[chromosome].append(line)
        except KeyError:
            kataegis_chromosomes[chromosome] = [line]
    return [kataegis_chromosomes, median(median_list)]


def process_kataegis_chromosomes(kataegis_chromosomes, sample_median):
    kataegis_hotspots = {}
    for key in kataegis_chromosomes.keys():
        kataegis_hotspots = find_hotspots(kataegis_chromosomes[key], sample_median, kataegis_hotspots, key)
    return kataegis_hotspots


def find_hotspots(kataegis_chromosome, sample_median, kataegis_hotspots, key):
    #if sample_median > 6000:
    #    sample_median = 6000
    sample_median = 6000
    kataegis_hotspots[key] = scan_window_by_six(kataegis_chromosome, sample_median)
    return kataegis_hotspots


def scan_window_by_six(kataegis_chromosome, sample_median):
    distances = find_intermutation_distance(kataegis_chromosome)
    found_spots = []
    i = 6
    while i < len(distances):
        if scan_list(distances[i-6: i], sample_median):
            found_spots.append(i)
            while scan_list(distances[i-6: i], sample_median):
                i += 1
        i += 1
    return found_spots


def find_intermutation_distance(kataegis_chromosome):
    distances = []
    for el in kataegis_chromosome:
        this_value = int(el.split("\t")[2])
        distances.append(this_value)
    return distances


def scan_list(list_to_scan, size):
    return sum(list_to_scan) < size


def median(lst):
    sorted_list = sorted(lst)
    index = (len(lst) - 1) // 2
    if len(lst) > 0:
        if len(lst) % 2:
            return sorted_list[index]
        else:
            return (sorted_list[index] + sorted_list[index + 1])/2.0
    else:
        return []