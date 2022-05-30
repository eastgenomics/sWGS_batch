# sWGS_batch
Collection of scripts to start batch jobs for sWGS.

The CNV caller used is WisecondorX: [https://academic.oup.com/nar/article/47/4/1605/5253050]

## How to run

The normal file is an optional file containing the ids of the samples which are the normal samples. It is a required input for running the full workflow. If not provided when running just the reference step, the script will assume that the files in the npz folder are the normal samples.

The sex file is a 2 column file containing the sample ids and their sex.

```bash
# align fastqs
python run.py align ${fastq_folder}

# downsample data
python run.py downsampling -ff ${flagstat_folder} -pf ${picard_folder} -c ${coverage} -o ${output_folder}

# cnv calling workflow
python run.py workflow -d ${downsampled_bam_folder} -n ${normal_file} -s ${sex_file} -b_npz ${binsize_npz} -b_ref ${binsize_ref} -o ${output_folder}

# npz conversion step
python run.py npz ${downsampled_bam_folder} -b ${binsize} -o ${output_folder}

# reference creation step
python run.py ref -f ${npz_folder} -b ${binsize} -n ${normal_file} -o ${output_folder}

# cnv calling step
python run.py -npz ${npz_folder} -r ${npz_reference} -s ${sex_file} -o ${output_folder}
```

## Output

The folder structure is similar to Dias:
- output
  - dias_single
  - downsampled folder (output folder is specified)
  - wisecondor output folder (output folder is specified)
    - npzs
    - output
      - plots folder
      - aberration bed
      - bins bed
      - segments bed
      - statistics txt
    - ref
