# xenium-filter-transcripts

This tool prepares the Xenium transcripts csv table to a format, so that it can be used with Baysor. This tools is adapted from the python script provided on the 10x [website](https://www.10xgenomics.com/cn/resources/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation).


## Preprocessing Xenium’s transcripts.csv.gz

For Baysor, the most relevant Xenium output file is transcripts.csv.gz. The file format is documented online [here](https://www.10xgenomics.com/cn/support/in-situ-gene-expression/documentation/steps/onboard-analysis/understanding-xenium-outputs#transcript-file).
Here are a few quick notes about some of its columns:

- **cell_id**: If a transcript is NOT associated with a cell, it will get -1. If it is associated with a cell, this column will have a positive integer cell ID.
- **x_location, y_location, z_location**: These three columns store the (x, y, z) coordinates of the transcript in microns. The origin is found in the upper-left corner of the bottom z-slice, as the following image illustrates.
- **feature_name**: The gene associated with the transcript.
- **qv**: Phred-scaled quality score of all decoded transcripts. It is up to the customer to determine a suitable Q-Score threshold when analyzing Xenium data. At 10x Genomics, we typically filter out transcripts with Q-Score < 20.

To preprocess Xenium’s transcript CSV file for Baysor, the script will perform the following steps:

- Filter out negative control transcripts.
- Filter out transcripts whose Q-Score falls below a specified threshold (default: 20).
- Baysor required that the **cell_id** is an integer. **UNASSIGNED** cells will be substituted with **0** and the encoded cell id will be reversed to an integer. The conversion is described [here](https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/xenium-outputs-zarr#cellID).

## Usage
   
    Usage: xenium-filter-transcripts [OPTIONS] <IN_FILE>
    
    Arguments:
    <IN_FILE>  The path to the transcripts.csv file produced by Xenium
    
    Options:
          --min-qv <MIN_QV>    The minimum Q-Score to pass filtering [default: 20]
          --min-x <MIN_X>      Only keep transcripts whose x-coordinate is greater than specified limit [default: 0]
          --max-x <MAX_X>      Only keep transcripts whose x-coordinate is less than specified limit.
                                 Xenium slide is <24000 microns in x and y. [default: 24000]
          --min-y <MIN_Y>      Only keep transcripts whose y-coordinate is greater than specified limit [default: 0]
          --max-y <MAX_Y>      Only keep transcripts whose y-coordinate is less than specified limit.
                                 Xenium slide is <24000 microns in x and y. [default: 24000]
          --nucleus-only       Only keep transcripts that are in the nucelus. All other transcripts will not be assigned to a cell
          --out-dir <OUT_DIR>  Path for outout file. The output file will be names after the following schema:  
                                 X{x-min}-{x-max}_Y{y-min}-{y-max}_filtered_transcripts_nucleus_only_{nucleus_only}.csv  
                                 E.g.: X0-24000_Y0-24000_filtered_transcripts_nucleus_only_false.csv [default: .]
      -h, --help               Print help
      -V, --version            Print version
  
