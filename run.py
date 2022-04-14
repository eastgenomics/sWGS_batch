import argparse

import downsample
import single
import utils


sWGS_project = "project-G8KJ7x84q42xKQJ95Fx2P7ZJ"
dias_single = "workflow-G8vy0v84q42qqBBf6GfgKZ89"

stages = {
    "sentieon_R1": "stage-Fy6fpk040vZZPPbq96Jb2KfK.reads_fastqgzs",
    "sentieon_R2": "stage-Fy6fpk040vZZPPbq96Jb2KfK.reads2_fastqgzs",
    "sentieon_sample": "stage-Fy6fpk040vZZPPbq96Jb2KfK.sample",
    "fastqc": "stage-Fy6fpV840vZZ0v6J8qBQYqZF.fastqs",
}


def main(args):
    cmd = args.which

    if cmd == "align":
        fastq_dict = single.make_fq_dict(args.fastq_folder)
        workflow = utils.setup_workflow(dias_single, sWGS_project)
        workflow_stage_info = utils.get_workflow_stage_info(workflow.id)
        stage_folders = utils.get_stage_output_folders(workflow_stage_info)
        single.run_workflow(workflow, fastq_dict, stages, stage_folders)

    elif cmd == "downsampling":
        app_handler = downsample.setup_downsampling_app()
        sample_mapping_numbers = downsample.get_mapping_numbers(
            args.flagstat_folder
        )
        sample_average_coverage = downsample.get_average_coverage(
            args.picard_folder
        )
        sample_downsampling_fraction = downsample.get_downsampling_fraction(
            sample_mapping_numbers, sample_average_coverage, args.coverage
        )
        sample_bams = downsample.gather_bams(sample_downsampling_fraction)
        downsample.start_downsampling_jobs(
            app_handler, sample_bams, sample_downsampling_fraction, args.output
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    align = subparsers.add_parser("align")
    align.add_argument(
        "fastq_folder", help="DNAnexus folder containing fastqs"
    )
    align.set_defaults(which="align")

    downsampling = subparsers.add_parser("downsampling")
    downsampling.add_argument(
        "flagstat_folder", help="Flagstat DNAnexus folder"
    )
    downsampling.add_argument("picard_folder", help="Picard DNAnexus folder")
    downsampling.add_argument(
        "-c", "--coverage", default=1, type=float,
        help="Wanted average coverage. Default=1X"
    )
    downsampling.add_argument(
        "-o", "--output", default=None,
        help="Output folder for the downsampling jobs"
    )
    downsampling.set_defaults(which="downsampling")

    args = parser.parse_args()
    main(args)
