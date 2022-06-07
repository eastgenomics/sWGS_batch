from collections import defaultdict
from pathlib import Path

import dxpy

import utils

""" setup function """


def setup_cnv_calling_app():
    """ Return wisecondorX app object for cnv calling

    Returns:
        DXApp: DXApp for wisecondorX
    """

    app_handler = dxpy.DXApp(name="eggd_wisecondorX")
    return app_handler


""" gathering functions """


def get_bams_and_bais(bam_folders):
    """ Get bams and bais given a list of folders

    Args:
        bam_folders (list): List of folders to look into for bams and bais

    Raises:
        Exception: if a bai is found but not the bam

    Returns:
        dict: Dict of sample ids and their bam/bai
    """

    print("Gathering bams and bais...")

    bams_bais = defaultdict(lambda: defaultdict(lambda: str))

    for bam_folder in bam_folders:
        bam_files = dxpy.find_data_objects(
            classname="file", name="*.bam", folder=bam_folder,
            name_mode="glob"
        )

        for bam in bam_files:
            bam = dxpy.DXFile(bam["id"])
            sample_id = "_".join(bam.name.split("_")[0:2])
            bams_bais[sample_id]["bam"] = [bam.id]

    for bam_folder in bam_folders:
        bai_files = dxpy.find_data_objects(
            classname="file", name="*.bai", folder=bam_folder,
            name_mode="glob"
        )

        for bai in bai_files:
            bai = dxpy.DXFile(bai["id"])
            sample_id = "_".join(bai.name.split("_")[0:2])

            if sample_id in bams_bais:
                bams_bais[sample_id]["bai"] = [bai.id]
            else:
                raise Exception(f"Found bai for {sample_id} but not a bam")

    assert bams_bais, "No bams or bais found"

    return bams_bais


def get_npzs_from_folders(npz_folders, npz_samples=[]):
    """ Get npz files from list of folders

    Args:
        npz_folders (list): List of folders to look into for npz
        npz_samples (list): List of samples to gather in folder

    Returns:
        dict: Dict of sample ids and their dnanexus links
    """

    print("Gathering npzs...")

    npzs = {}
    filtered_npzs = {}

    for npz_folder in npz_folders:
        npz_files = dxpy.find_data_objects(
            classname="file", name="*.npz", folder=npz_folder,
            name_mode="glob"
        )

        for npz in npz_files:
            npz_file = dxpy.DXFile(npz["id"])
            npz_file_name = Path(npz_file.name).stem
            sample_id = "_".join(npz_file_name.split("_")[0:2])
            npzs[sample_id] = [npz_file.id]

    if npz_samples:
        for sample_id in npzs:
            if sample_id in npz_samples:
                filtered_npzs[sample_id] = npzs[sample_id]
    else:
        filtered_npzs = npzs

    assert filtered_npzs, "No npz files found"

    return filtered_npzs


def get_normal_samples(npzs, normal_samples):
    """Given a dict of sample ids and their npz file and a list of sample ids,
    create a new dict with only the samples in the list

    Args:
        npzs (dict): Dict of sample ids and their npz file
        normal_samples (list): List of sample ids

    Returns:
        dict: Dict of sample ids and their npz file
    """

    normals = {}

    for sample_id in npzs:
        if sample_id in normal_samples:
            normals[sample_id] = npzs[sample_id]

    return normals


def get_sex_of_sample(sample_id, sample_sexes):
    """ Given a sample id and dict of sample ids and their sex, return the sex.
    If the sample is not found in the dict, return None

    Args:
        sample_id (str): Sample id
        sample_sexes (dict): Dict of sample ids and their sex

    Returns:
        str: One letter string to define sex
    """

    if sample_id in sample_sexes:
        return sample_sexes[sample_id]

    return None


""" parsing functions """


def parse_sex_file(sex_file):
    """ Parse file containing sample ids and the sex of those samples

    Args:
        sex_file (str): Path to the sex file

    Returns:
        dict: Dict of sample ids and their sex
    """

    sample_sexes = {}

    with open(sex_file) as f:
        for line in f:
            sample, sex = line.strip().split()
            sample_sexes[sample] = sex

    return sample_sexes


def parse_normal_sample_file(file):
    """ Parse file containing normal sample ids

    Args:
        file (str): Path to normal sample file

    Returns:
        list: List of normal samples
    """

    normal_samples = []

    with open(file) as f:
        for line in f:
            normal_samples.append(line.strip())

    return normal_samples


""" job functions """


def convert_npz(
    app_handler, bam_folder, binsize, out_folder, out_npzs_folder=None,
    ref_binsize=None
):
    """ Start jobs for converting bams to npz

    Args:
        app_handler (DXApp): DXApp object for wisecondorX
        bam_folder (list): List of folders to look for bams and bais into
        binsize (int): Binsize for npz
        out_folder (str): Out folder for DNAnexus
        out_npzs_folder (str, optional): Out folder for the npzs specifically.
                                Defaults to None
        ref_binsize (int, optional): Binsize for reference. Defaults to None

    Returns:
        dict: Dict of sample ids and their job objects
    """

    if out_npzs_folder:
        folder = out_npzs_folder
    else:
        folder = f"{out_folder}/npzs"

    sample_bams_bais = get_bams_and_bais(bam_folder)

    jobs = {}

    print("Setting up conversion jobs...")

    for sample_id in sample_bams_bais:
        if ref_binsize:
            job_name = (
                f"{app_handler.name} - Workflow {binsize}kb/{ref_binsize}kb - "
                f"Convert npz ({binsize}) - {sample_id}"
            )
        else:
            job_name = (
                f"{app_handler.name} - Convert npz ({binsize}) - {sample_id}"
            )

        inputs = {}
        inputs["bam"] = utils.create_dnanexus_links(
            sample_bams_bais[sample_id]["bam"], app_handler, "bam"
        )
        inputs["bai"] = utils.create_dnanexus_links(
            sample_bams_bais[sample_id]["bai"], app_handler, "bai"
        )
        inputs["binsize_convert"] = binsize
        inputs["convert_npz"] = True
        jobs[sample_id] = app_handler.run(inputs, folder=folder, name=job_name)

    print("Conversion jobs started...")

    return jobs


def create_ref(
    app_handler, binsize, out_folder, normal_file=None, npz_folders=None,
    npz_jobs=None, npz_binsize=None
):
    """ Start reference creation job

    Args:
        app_handler (DXApp): DXApp for wisecondorX
        binsize (int): Binsize for reference
        out_folder (str): Out folder in DNAnexus
        normal_file (str, optional): Path to file containing normal sample ids.
                                    Defaults to None.
        npz_folders (list, optional): Folders containing the npz files.
                                    Defaults to None.
        npz_jobs (dict, optional): Dict containing sample ids and their
                                    conversion job. Defaults to None.
        npz_binsize (int, optional): Npz binsize. Defaults to None

    Raises:
        Exception: check if npz folder or npz job is passed

    Returns:
        DXJob: DXJob for the reference creation job
    """

    npzs = {}

    # prepare the npzs for the ref job
    if npz_folders:
        npzs = get_npzs_from_folders(npz_folders)
        npzs = {
            sample_id: utils.create_dnanexus_links(
                npzs[sample_id], app_handler, "npz"
            )
            for sample_id in npzs
        }
    elif npz_jobs:
        for sample_id in npz_jobs:
            npzs[sample_id] = npz_jobs[sample_id].get_output_ref(
                "wisecondorx_output"
            )
    else:
        raise Exception("No npz folder or npz jobs passed")

    # filter out to get only the normal samples
    if normal_file:
        normal_samples = parse_normal_sample_file(normal_file)
        input_normal_samples = get_normal_samples(npzs, normal_samples)
    else:
        # if no normal file is passed, assume that the folder given contains
        # only normals
        input_normal_samples = [npz.id for npz in npzs]

    inputs = {}
    inputs["npz"] = list(input_normal_samples.values())
    inputs["binsize_newref"] = binsize
    inputs["create_ref"] = True

    if npz_binsize:
        job_name = (
            f"{app_handler.name} - Workflow {npz_binsize}kb/{binsize}kb - "
            f"Create ref ({binsize})"
        )
    else:
        job_name = f"{app_handler.name} - Create ref ({binsize})"

    assert inputs["npz"], "No normal npz files found."

    print("Setting up reference job...")

    # setup appropriate instance type
    if binsize > 2000 and binsize <= 5000:
        instance_type = "mem1_ssd1_v2_x36"
    elif binsize > 1000 and binsize <= 2000:
        instance_type = "mem1_ssd1_v2_x72"
    elif binsize <= 1000:
        instance_type = "mem2_ssd1_v2_x64"
    else:
        instance_type = "mem1_ssd1_v2_x8"

    if npz_jobs:
        ref_job = app_handler.run(
            inputs, folder=f"{out_folder}/ref", name=job_name,
            depends_on=list(npz_jobs.values()),
            instance_type=instance_type
        )
    else:
        ref_job = app_handler.run(
            inputs, folder=f"{out_folder}/ref",
            instance_type=instance_type, name=job_name,
        )

    print("Reference job started...")

    return ref_job


def call_cnvs(
    app_handler, sex_file, out_folder, npz_jobs=None, ref_job=None,
    npz_samples=[], npz_folders=None, ref_path=None, npz_binsize=None,
    ref_binsize=None
):
    """ Start cnv calling jobs

    Args:
        app_handler (DXApp): DXApp for wisecondorX
        sex_file (str): Path to file containing sample ids and their sex
        out_folder (str): Out folder in DNAnexus
        npz_jobs (dict, optional): Dict of sample ids and their npz jobs.
                                    Defaults to None.
        ref_job (DXJob, optional): DXJob for the creation of reference.
                                    Defaults to None.
        npz_samples (str, optional): DNAnexus path to find npz files in.
                                    Defaults to None.
        npz_folders (str, optional): DNAnexus path to find npz files in.
                                    Defaults to None.
        ref_path (str, optional): DNAnexus path to reference file.
                                    Defaults to None.
        npz_binsize (int, optional): Npz binsize. Defaults to None
        ref_binsize (int, optional): Ref binsize. Defaults to None

    Raises:
        Exception: Check if npz folder or npz jobs passed
        Exception: Check if ref folder or ref job passed
    """

    sample_sexes = parse_sex_file(sex_file)

    npzs = {}

    # prepare npzs for the cnv calling jobs
    if npz_folders:
        npzs = get_npzs_from_folders(npz_folders, npz_samples)
        npzs = {
            sample_id: utils.create_dnanexus_links(
                npzs[sample_id], app_handler, "npz"
            )
            for sample_id in npzs
        }
    elif npz_jobs:
        for sample_id in npz_jobs:
            npzs[sample_id] = npz_jobs[sample_id].get_output_ref(
                "wisecondorx_output"
            )
    else:
        raise Exception("No npz folder or npz jobs passed")

    # prepare the reference for the cnv calling jobs
    if ref_path:
        assert ref_path.startswith("/"), (
            "DNAnexus paths need to start with a /"
        )
        folder_path = str(Path(ref_path).parent)
        ref_name = str(Path(ref_path).name)
        ref_object = dxpy.find_one_data_object(
            folder=folder_path, name=ref_name
        )
        ref = utils.create_dnanexus_links(
            [ref_object["id"]], app_handler, "reference"
        )
    elif ref_job:
        ref = ref_job.get_output_ref("wisecondorx_output")
    else:
        raise Exception("No ref folder or ref job passed")

    print("Setting up cnv calling jobs...")

    for sample_id in npzs:
        inputs = {}
        inputs["npz"] = npzs[sample_id]
        inputs["reference"] = ref
        inputs["variant_calling"] = True

        sex = get_sex_of_sample(sample_id, sample_sexes)

        if sex:
            inputs["sex"] = sex

        if npz_binsize and ref_binsize:
            job_name = (
                f"{app_handler.name} - Workflow {npz_binsize}kb/{ref_binsize}"
                f"kb - CNV calling ({sex}) - {sample_id}"
            )
        else:
            job_name = (
                f"{app_handler.name} - CNV calling ({sex}) - {sample_id}"
            )

        if npz_jobs and ref_job:
            jobs = [npz_jobs[sample_id]] + [ref_job]
            app_handler.run(
                inputs, folder=f"{out_folder}/output", depends_on=jobs,
                name=job_name
            )
        else:
            app_handler.run(
                inputs, folder=f"{out_folder}/output", name=job_name
            )

    print("Cnv calling jobs started...")


""" workflow function """


def run_cnv_calling(
    app_handler, bam_folder, normal_file, sex_file, npz_binsize, ref_binsize,
    out_folder, out_npz_folder
):
    """ Run the full "workflow" for cnv calling using wisecondorX

    Args:
        app_handler (DXApp): DXApp for wisecondorX
        bam_folder (list): List of bam folders to look into for bams and bais
        normal_file (str): Path to normal file
        sex_file (str): Path to sex file
        npz_binsize (int): Binsize for convert npz jobs
        ref_binsize (int): Binsize for create reference job
        out_folder (str): DNAnexus out folder
        out_npz_folder (str): Npz DNAnexus out folder
    """

    npz_jobs = convert_npz(
        app_handler, bam_folder, npz_binsize, out_folder, out_npz_folder,
        ref_binsize=ref_binsize
    )
    ref_job = create_ref(
        app_handler, ref_binsize, out_folder, normal_file=normal_file,
        npz_jobs=npz_jobs, npz_binsize=npz_binsize
    )
    call_cnvs(
        app_handler, sex_file, out_folder, npz_jobs=npz_jobs, ref_job=ref_job,
        npz_binsize=npz_binsize, ref_binsize=ref_binsize
    )


def run_from_npzs(
    app_handler, npz_folders, binsize, normal_file, sex_file, out_folder
):
    """ Run the cnv calling from already generated npzs

    Args:
        app_handler (DXApp): App handler for WisecondorX
        npz_folders (list): List of npz folders
        binsize (int): Binsize for reference creation
        normal_file (str): Path to normal file
        sex_file (str): Path to sex file
        out_folder (str): DNAnexus out folder
    """

    ref_job = create_ref(
        app_handler, binsize, out_folder, normal_file=normal_file,
        npz_folders=npz_folders
    )
    call_cnvs(
        app_handler, sex_file, out_folder, npz_folders=npz_folders,
        ref_job=ref_job
    )
