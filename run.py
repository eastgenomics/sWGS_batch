import argparse

import downsample
import single
import utils
import cnv_calling


sWGS_project = "project-G8KJ7x84q42xKQJ95Fx2P7ZJ"
dias_single = "workflow-G8vy0v84q42qqBBf6GfgKZ89"

stages = {
    "sentieon_R1": "stage-Fy6fpk040vZZPPbq96Jb2KfK.reads_fastqgzs",
    "sentieon_R2": "stage-Fy6fpk040vZZPPbq96Jb2KfK.reads2_fastqgzs",
    "sentieon_sample": "stage-Fy6fpk040vZZPPbq96Jb2KfK.sample",
    "fastqc": "stage-Fy6fpV840vZZ0v6J8qBQYqZF.fastqs",
}


def main(args):
    cmd = args.cmd

    if cmd == "align":
        fastq_dict = single.make_fq_dict(args.fastq_folder)
        workflow = utils.setup_workflow(dias_single, sWGS_project)
        workflow_stage_info = utils.get_workflow_stage_info(workflow.id)
        stage_folders = utils.get_stage_output_folders(workflow_stage_info)
        single.run_workflow(workflow, fastq_dict, stages, stage_folders)

    elif cmd == "downsampling":
        picard_app_handler = downsample.setup_downsampling_app()
        samtools_app_handler = downsample.setup_indexing_app()
        sample_mapping_numbers = downsample.get_mapping_numbers(
            args.flagstat_folders
        )
        sample_average_coverage = downsample.get_average_coverage(
            args.picard_folders
        )
        sample_downsampling_fraction = downsample.get_downsampling_fraction(
            sample_mapping_numbers, sample_average_coverage, args.coverage
        )
        sample_bams = downsample.gather_bams(sample_downsampling_fraction)
        out_bams = downsample.start_downsampling_jobs(
            picard_app_handler, samtools_app_handler, sample_bams,
            sample_downsampling_fraction, args.output
        )

    elif cmd == "cnv_calling":
        cnv_cmd = args.cnv_op
        app_handler = cnv_calling.setup_cnv_calling_app()

        if cnv_cmd == "workflow":
            cnv_calling.run_cnv_calling(
                app_handler, args.downsampled_bam_folder, args.normal_file,
                args.sex_file, args.binsize_npz, args.binsize_ref,
                args.out_folder
            )

        if cnv_cmd == "npz":
            cnv_calling.convert_npz(
                app_handler, args.downsampled_bam_folder, args.binsize,
                args.out_folder
            )

        if cnv_cmd == "ref":
            cnv_calling.create_ref(
                app_handler, args.binsize, args.out_folder,
                normal_file=args.normal_file, npz_folder=args.npz_folder
            )

        if cnv_cmd == "cnv":
            cnv_calling.call_cnvs(
                app_handler, args.sex_file, args.out_folder, 
                npz_folder=args.npz_folder, ref_path=args.npz_reference
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    align = subparsers.add_parser("align")
    align.add_argument(
        "fastq_folder", nargs="+", help="DNAnexus folder containing fastqs"
    )
    align.set_defaults(which="align")

    downsampling = subparsers.add_parser("downsampling")
    downsampling.add_argument(
        "-ff", "--flagstat_folders", nargs="+", help="Flagstat DNAnexus folder"
    )
    downsampling.add_argument(
        "-pf", "--picard_folders", nargs="+", help="Picard DNAnexus folder"
    )
    downsampling.add_argument(
        "-c", "--coverage", default=1, type=float,
        help="Wanted average coverage. Default=1X"
    )
    downsampling.add_argument(
        "-o", "--output", default=None,
        help="Output folder for the downsampling jobs"
    )
    downsampling.set_defaults(which="downsampling")

    cnv_calling_parser = subparsers.add_parser("cnv_calling")
    cnv_subparser = cnv_calling_parser.add_subparsers(dest="cnv_op")

    full_workflow = cnv_subparser.add_parser(
        "workflow",
        help="Run the cnv calling without having to run individual commands",
        description=(
            "Run the cnv calling without having to run individual commands"
        )
    )
    full_workflow.add_argument(
        "-d", "--downsampled_bam_folder", nargs="+",
        help="Downsampled bam DNAnexus folder"
    )
    full_workflow.add_argument(
        "-n", "--normal_file", default=None,
        help="File containing normal samples"
    )
    full_workflow.add_argument(
        "-s", "--sex_file", help="File containing the sex of the samples"
    )
    full_workflow.add_argument(
        "-b_npz", "--binsize_npz", type=int,
        help="Binsize for converting bams into npz"
    )
    full_workflow.add_argument(
        "-b_ref", "--binsize_ref", type=int, help="Binsize for reference"
    )
    full_workflow.add_argument(
        "-o", "--out_folder",
        help=(
            "Output folder where output of the apps will dump their own out "
            "folders"
        )
    )

    npz_generation = cnv_subparser.add_parser(
        "npz", help="Individual step to generate npz files from BAM and BAI",
        description="Individual step to generate npz files from BAM and BAI"
    )
    npz_generation.add_argument(
        "downsampled_bam_folder", nargs="+",
        help="Downsampled bam DNAnexus folder"
    )
    npz_generation.add_argument("-b", "--binsize", help="Binsize for npz")
    npz_generation.add_argument(
        "-o", "--out_folder",
        help=(
            "Output folder where output of the app will dump its own out "
            "folder"
        )
    )

    ref_generation = cnv_subparser.add_parser(
        "ref", help="Individual step to generate npz reference file",
        description="Individual step to generate npz reference file"
    )
    ref_generation.add_argument(
        "-f", "--npz_folder", nargs="+",
        help="Folder containing npz files in DNAnexus"
    )
    ref_generation.add_argument(
        "-b", "--binsize", type=int, help="Binsize for the reference file"
    )
    ref_generation.add_argument(
        "-n", "--normal_file", default=None,
        help="File containing normal samples"
    )
    ref_generation.add_argument(
        "-o", "--out_folder",
        help=(
            "Output folder where output of the app will dump its own out "
            "folder"
        )
    )

    cnv_generation = cnv_subparser.add_parser(
        "cnv", help="Individual step to generate cnv calling output",
        description="Individual step to generate cnv calling output"
    )
    cnv_generation.add_argument(
        "-npz", "--npz_folder", nargs="+", help="Folder containing npz files"
    )
    cnv_generation.add_argument(
        "-r", "--npz_reference", help="DNAnexus path to npz reference"
    )
    cnv_generation.add_argument(
        "-s", "--sex_file", help="File containing the sex of the samples"
    )
    cnv_generation.add_argument(
        "-o", "--out_folder",
        help=(
            "Output folder where output of the app will dump its own out "
            "folder"
        )
    )

    args = parser.parse_args()
    main(args)
