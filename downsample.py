import dxpy

import utils


def setup_downsampling_app():
    """ Setup DXApp Picard downsample

    Returns:
        dxpy.DXApp: DXApp object
    """

    app_handler = dxpy.DXApp(name="picard_downsample")
    return app_handler


def get_mapping_numbers(flagstat_folder):
    """ Get the number of reads mapped for each sample

    Args:
        flagstat_folder (str): DNAnexus folder where flagstats files are

    Returns:
        dict: Dict of sample name with number of mapped reads
    """

    flagstats = {}

    flagstat_files = dxpy.find_data_objects(
        classname="file", name="*.flagstat", folder=flagstat_folder,
        name_mode="glob"
    )

    for flagstat_file in flagstat_files:
        flagstat_file_name = dxpy.DXFile(flagstat_file["id"]).name
        sample_name = "_".join(flagstat_file_name.split("_")[:2])

        with dxpy.open_dxfile(flagstat_file["id"]) as f:
            for line in f:
                line = line.strip().split()

                if "mapped" in line:
                    nb_mapped = line[0]
                    flagstats[sample_name] = nb_mapped
                    break

    return flagstats


def get_average_coverage(picard_folder):
    """ Get mean coverage number obtained through picard WGS QC

    Args:
        picard_folder (str): Folder where the Picard QC is located

    Returns:
        dict: Dict of sample name to average coverage
    """

    picard = {}

    picard_files = dxpy.find_data_objects(
        classname="file", name="*.wgs_stats.tsv", folder=picard_folder,
        name_mode="glob"
    )

    for picard_file in picard_files:
        picard_file_name = dxpy.DXFile(picard_file["id"]).name
        sample_name = "_".join(picard_file_name.split("_")[:2])

        data_bool = False

        with dxpy.open_dxfile(picard_file["id"]) as f:
            for line in f:
                if "GENOME_TERRITORY" in line:
                    data_bool = True
                    headers = line.strip().split()
                    continue

                if data_bool:
                    data = line.strip().split()
                    break

        for header, data_ele in zip(headers, data):
            if header == "MEAN_COVERAGE":
                picard[sample_name] = data_ele

    return picard


def get_downsampling_fraction(flagstat_data, picard_data, coverage):
    """ Calculate downsampling fraction parameter to be passed to Picard
        Downsample

    Args:
        flagstat_data (dict): Dict from sample to number of mapped reads
        picard_data (dict): Dict from sample to average coverage
        coverage (float): Desired coverage for final downsampled BAM

    Returns:
        dict: Dict from sample to downsampling fraction
    """

    downsample_fraction = {}

    for sample in flagstat_data:
        nb_reads_for_wanted_coverage = float(flagstat_data[sample]) / float(picard_data[sample]) * float(coverage)
        downsampling_fraction = float(nb_reads_for_wanted_coverage) / float(flagstat_data[sample])
        downsample_fraction[sample] = downsampling_fraction

    return downsample_fraction


def gather_bams(sample_dict):
    """ Get BAM file ids for samples in the project

    Args:
        sample_dict (dict): Dict from sample to downsampling fraction

    Returns:
        dict: Dict from sample to BAM file ids
    """

    bam_dict = {}

    for sample in sample_dict:
        bam_file = dxpy.find_one_data_object(
            more_ok=False, classname="file", name=f"{sample}_*.bam",
            name_mode="glob"
        )

        bam_dict[sample] = [bam_file["id"]]

    return bam_dict


def start_downsampling_jobs(
    app_handler, bam_dict, downsampling_dict, out_folder
):
    """ For every sample, start a downsampling job

    Args:
        app_handler (dxpy.DXApp): Picard downsample DXApp
        bam_dict (dict): Dict containing samples and their BAM file ids
        downsampling_dict (dict): Dict containing samples and the downsampling fraction
        out_folder (str): Output folder in DNAnexus
    """

    if out_folder is None:
        out_folder = f"{app_handler.describe()['name']}_v1.1.0"

    inputs = {}

    for sample in bam_dict:
        inputs["sorted_bam"] = utils.create_dnanexus_links(bam_dict[sample])
        inputs["fraction"] = downsampling_dict[sample]
        app_handler.run(inputs, tags=["1.1.0"], folder=out_folder)
