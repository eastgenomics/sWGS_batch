import os

import dxpy

import utils


def setup_workflow():
    """ Setup single workflow for sWGS

    Returns:
        DXWorkflow: DXWorkflow object for single workflow
    """

    sWGS_project = "project-G8KJ7x84q42xKQJ95Fx2P7ZJ"
    dias_single = "workflow-G8vy0v84q42qqBBf6GfgKZ89"
    workflow = dxpy.DXWorkflow(dias_single, project=sWGS_project)
    return workflow


def make_fq_dict(path):
    """ Match fastq files located in the given DNAnexus path

    Args:
        path (str): DNAnexus path

    Returns:
        dict: Dict containing R1 and R2 fastqs
    """

    # get all the files ending by fastq.gz in the DNAnexus path

    fastq_dict = {}

    print("Gathering fastqs...")

    try:
        fastqs = dxpy.find_data_objects(
            name="*fastq.gz", folder=path, name_mode="glob"
        )
        next(fastqs)
    except StopIteration as e:
        print(
            f"No fastqs found in {os.environ['DX_PROJECT_CONTEXT_ID']}{path}"
        )
    else:
        fastqs = dxpy.find_data_objects(
            name="*fastq.gz", folder=path, name_mode="glob"
        )

    for fastq in fastqs:
        fastq_obj = dxpy.DXFile(fastq["id"])

        if not fastq_obj.name.startswith("Undetermined"):
            sample_id = "_".join(fastq_obj.name.split("_")[:2])

            if sample_id == "Undetermined":
                continue

            read_num = None
            if "_R1_" in fastq_obj.name:
                read_num = "R1"
            elif "_R2_" in fastq_obj.name:
                read_num = "R2"

            assert read_num, (
                "Unable to determine read number (R1 or R2) for fastq "
                f"{fastq_obj.name}"
            )

            # Make a new dict entry for sample if not present
            fastq_dict.setdefault(sample_id, {
                "R1": [],
                "R2": []}
            )

            # Add fastq filename and file_id to appropriate place in dict.
            # We add both id and name because we need to sort by name later
            fastq_dict[sample_id].setdefault(read_num, []).append(
                (fastq_obj.name, fastq_obj.id)
            )

    # Sort fastq lists so that the fastq at pos n in R1 list
    # is paired with the fastq at pos n in R2 list
    # Once the sort is complete we remove the filename from the dict
    # since it was only there to enable the sort
    for sample in fastq_dict:
        for read in ["R1", "R2"]:
            # sort tuple on first element i.e. filename
            # retain only file id i.e. second element
            sorted_fastq_list = [
                x[1]
                for x in sorted(fastq_dict[sample][read], key=lambda x: x[0])
            ]
            fastq_dict[sample][read] = sorted_fastq_list

    return fastq_dict


def setup_inputs_workflow(stages, sample, reads):
    """ Setup the inputs for the single workflow

    Args:
        stages (dict): Dict of the stages and their stage id input
        sample (str): Sample id
        reads (dict): Dict of reads

    Returns:
        dict: Dict of stage ids and their dnanexus links inputs
    """

    inputs = {}

    inputs[stages["sentieon_R1"]] = utils.create_dnanexus_links(reads["R1"])
    inputs[stages["sentieon_R2"]] = utils.create_dnanexus_links(reads["R2"])
    inputs[stages["fastqc"]] = utils.create_dnanexus_links(
        reads["R1"] + reads["R2"]
    )
    inputs[stages["sentieon_sample"]] = sample

    return inputs


def run_workflow(workflow, fastq_dict, stages, stage_folders):
    """ Start the single workflow jobs

    Args:
        workflow (DXWorkflow): Workflow object
        fastq_dict (dict): Dict of sample and their fastqs
        stages (dict): Dict of stages and their input name
        stage_folders (dict): Dict of stages and their output folder
    """

    run_datetime = utils.get_run_datetime()

    for sample, reads in fastq_dict.items():
        assert(len(reads["R1"]) == len(reads["R2"])), \
            "Mismatched number of R1/R2 fastqs for {}".format(sample)

        inputs = setup_inputs_workflow(stages, sample, reads)
        workflow.run(
            inputs, folder=f"/output_{run_datetime}",
            stage_folders=stage_folders
        )
        print(f"Started Sentieon for {sample}")
