from collections import defaultdict
import dxpy

import utils

""" setup function """


def setup_cnv_calling_app():
    app_handler = dxpy.DXApp(name="eggd_wisecondorX")
    return app_handler


""" gathering functions """


def get_bams_and_bais(bam_folders):
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


def get_npzs_from_folder(npz_folders):
    print("Gathering npzs...")

    npzs = {}

    for npz_folder in npz_folders:
        npz_files = dxpy.find_data_objects(
            classname="file", name="*.npz", folder=npz_folder,
            name_mode="glob"
        )

        for npz in npz_files:
            npz_file = dxpy.DXFile(npz["id"])
            sample_id = "_".join(npz_file.name.split("_")[0:2])
            dnanexus_link = utils.create_dnanexus_links(npz_file.id)
            npzs[sample_id] = dnanexus_link

    return npzs


def get_ref_path_from_job(ref_job):
    return f"{ref_job.describe()['folder']}/ref.npz"


def get_normal_samples(npzs, normal_samples):
    normals = {}

    for sample_id in npzs:
        if sample_id in normal_samples:
            normals[sample_id] = npzs[sample_id]

    return normals


def get_sex_of_sample(sample_id, sample_sexes):
    if sample_id in sample_sexes:
        return sample_sexes[sample_id]

    return None


def get_ref_output_from_job(job, output_field):
    return job.get_output_ref(output_field)


""" parsing functions """


def parse_sex_file(sex_file):
    sample_sexes = {}

    with open(sex_file) as f:
        for line in f:
            sample, sex = line.strip().split()
            sample_sexes[sample] = sex

    return sample_sexes


def parse_normal_sample_file(file):
    normal_samples = []

    with open(file) as f:
        for line in f:
            normal_samples.append(line.strip())

    return normal_samples


""" job functions """


def convert_npz(app_handler, bam_folder, binsize, out_folder):

    sample_bams_bais = get_bams_and_bais(bam_folder)

    jobs = {}
    inputs = {}

    print("Setting up conversion jobs...")

    for sample_id in sample_bams_bais:
        inputs["bam"] = [
            utils.create_dnanexus_links(sample_bams_bais[sample_id]["bam"])
        ]
        inputs["bai"] = [
            utils.create_dnanexus_links(sample_bams_bais[sample_id]["bai"])
        ]
        inputs["binsize_convert"] = binsize
        inputs["convert_npz"] = True
        job_name = (
            f"{app_handler.name} - Convert npz ({binsize}) - {sample_id}"
        )
        jobs[sample_id] = app_handler.run(
            inputs, folder=f"{out_folder}/npzs", name=job_name
        )

    print("Conversion jobs started...")

    return jobs


def create_ref(
    app_handler, binsize, out_folder, normal_file=None, npz_folder=None,
    npz_jobs=None,
):
    npzs = {}

    if npz_folder:
        npzs = get_npzs_from_folder(npz_folder)
    elif npz_jobs:
        for sample_id in npz_jobs:
            npzs[sample_id] = get_ref_output_from_job(
                npz_jobs[sample_id], "wisecondorx_output"
            )
    else:
        raise Exception("No npz folder or npz jobs passed")

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
    job_name = f"{app_handler.name} - Create ref ({binsize})"

    print("Setting up reference job...")

    if npz_jobs:
        ref_job = app_handler.run(
            inputs, folder=f"{out_folder}/ref", name=job_name,
            depends_on=list(npz_jobs.values())
        )
    else:
        ref_job = app_handler.run(inputs, folder=f"{out_folder}/ref")

    print("Reference job started...")

    return ref_job


def call_cnvs(
    app_handler, sex_file, out_folder, npz_jobs=None, ref_job=None,
    npz_folder=None, ref_path=None
):
    sample_sexes = parse_sex_file(sex_file)

    npzs = {}

    if npz_folder:
        npzs = get_npzs_from_folder(npz_folder)
    elif npz_jobs:
        for sample_id in npz_jobs:
            npzs[sample_id] = get_ref_output_from_job(
                npz_jobs[sample_id], "wisecondorx_output"
            )
    else:
        raise Exception("No npz folder or npz jobs passed")

    if ref_path:
        ref = utils.create_dnanexus_links(dxpy.DXFile(ref_path).id)
    elif ref_job:
        ref = get_ref_output_from_job(
            ref_job, "wisecondorx_output"
        )
        print(ref)
    else:
        raise Exception("No ref folder or ref job passed")

    inputs = {}

    print("Setting up cnv calling jobs...")

    for sample_id in npzs:
        inputs["npz"] = [npzs[sample_id]]
        inputs["reference"] = ref
        inputs["variant_calling"] = True

        sex = get_sex_of_sample(sample_id, sample_sexes)

        if sex:
            inputs["sex"] = sex

        job_name = f"{app_handler.name} - CNV calling ({sex}) - {sample_id}"

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
    out_folder
):
    npz_jobs = convert_npz(app_handler, bam_folder, npz_binsize, out_folder)
    ref_job = create_ref(
        app_handler, ref_binsize, out_folder, normal_file=normal_file,
        npz_jobs=npz_jobs
    )
    call_cnvs(
        app_handler, sex_file, out_folder, npz_jobs=npz_jobs, ref_job=ref_job
    )
