use clap::Parser;
use serde::{Deserialize, Serialize};
use std::{error::Error, fs::File};

/// Filter transcripts from transcripts.csv based on Q-Score threshold
/// and upper bounds on x and y coordinates. Remove negative controls.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// The path to the transcripts.csv file produced
    /// by Xenium.
    in_file: String,

    /// The minimum Q-Score to pass filtering.
    #[arg(long, default_value_t = 20.0)]
    min_qv: f32,

    /// Only keep transcripts whose x-coordinate is greater than specified limit.
    #[arg(long, default_value_t = 0.0)]
    min_x: f32,

    /// Only keep transcripts whose x-coordinate is less than specified limit.
    ///   Xenium slide is <24000 microns in x and y.
    #[arg(long, default_value_t = 24000.0, verbatim_doc_comment)]
    max_x: f32,

    /// Only keep transcripts whose y-coordinate is greater than specified limit.
    #[arg(long, default_value_t = 0.0)]
    min_y: f32,

    /// Only keep transcripts whose y-coordinate is less than specified limit.
    ///   Xenium slide is <24000 microns in x and y.
    #[arg(long, default_value_t = 24000.0, verbatim_doc_comment)]
    max_y: f32,

    /// Only keep transcripts that are in the nucelus.
    /// All other transcripts will not be assigned to a cell.
    #[arg(long, default_value_t = false)]
    nucleus_only: bool,

    /// Path for outout file. The output file will be names after the following schema:  
    ///   X{x-min}-{x-max}_Y{y-min}-{y-max}_filtered_transcripts_nucleus_only_{nucleus_only}.csv  
    ///   E.g.: X0-24000_Y0-24000_filtered_transcripts_nucleus_only_false.csv
    #[arg(long, default_value = ".", verbatim_doc_comment)]
    out_dir: String,
}

#[derive(Debug, Deserialize, Serialize)]
struct Transcript {
    transcript_id: usize,
    cell_id: String,
    overlaps_nucleus: usize,
    feature_name: String,
    x_location: f32,
    y_location: f32,
    z_location: f32,
    qv: f32,
    fov_name: String,
    nucleus_distance: f32,
}

struct CellID {
    cell_id_prefix: usize,
    dataset_suffix: usize,
}

impl Default for Transcript {
    fn default() -> Transcript {
        Transcript {
            transcript_id: 1,
            cell_id: "a".to_string(),
            overlaps_nucleus: 0,
            feature_name: "gene".to_string(),
            x_location: 100.0,
            y_location: 100.0,
            z_location: 100.0,
            qv: 25.0,
            fov_name: "FOV1".to_string(),
            nucleus_distance: 1.0,
        }
    }
}

impl Transcript {
    fn new() -> Self {
        Default::default()
    }
}

pub fn run(args: Args) -> Result<(), Box<dyn Error>> {
    let file = File::open(&args.in_file)?;
    let mut rdr = csv::Reader::from_reader(file);

    let mut wtr = create_out_file(&args)?;

    for result in rdr.deserialize() {
        let mut record: Transcript = result?;

        // Check the coordinates and qv values
        let keep = filter_transcripts(
            &record,
            args.min_x,
            args.max_x,
            args.min_y,
            args.max_y,
            args.min_qv,
        );
        if keep {
            // remove cell assigment if not in nucleous and
            // if arg nucleus_only
            if args.nucleus_only && record.overlaps_nucleus == 0 {
                record.cell_id = "0".to_string();
            }

            // Change cell_id of cell-free transcripts from UNASSIGNED to 0
            if record.cell_id == "UNASSIGNED" {
                record.cell_id = "0".to_string();
            }

            // Decode cell_ids
            if record.cell_id != "0" {
                record.cell_id = decode_cell_id(&record.cell_id).cell_id_prefix.to_string();
            }

            wtr.serialize(&record)?;
        }
    }
    wtr.flush()?;
    Ok(())
}

fn create_out_file(args: &Args) -> Result<csv::Writer<File>, Box<dyn Error>> {
    std::fs::create_dir_all(&args.out_dir)?;

    let outfile = format!(
        "{}/X{}-{}_Y{}-{}_filtered_transcripts_nucleus_only_{}.csv",
        args.out_dir, args.min_x, args.max_x, args.min_y, args.max_y, args.nucleus_only
    );

    let wtr = csv::Writer::from_path(outfile)?;
    Ok(wtr)
}

fn filter_transcripts(
    t: &Transcript,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    min_qc: f32,
) -> bool {
    let keep = t.x_location >= min_x
        && t.x_location <= max_x
        && t.y_location >= min_y
        && t.y_location <= max_y
        && t.qv >= min_qc
        && !exclude_feature_names(t);
    return keep;
}

fn exclude_feature_names(t: &Transcript) -> bool {
    return t.feature_name.starts_with("NegControlProbe_")
        || t.feature_name.starts_with("antisense_")
        || t.feature_name.starts_with("NegControlCodeword_")
        || t.feature_name.starts_with("BLANK_");
}

fn decode_cell_id(cell_id: &str) -> CellID {
    let mut parts = cell_id.split("-");
    let shifted_hex_digits = parts.next().unwrap();
    let dataset_suffix: usize = parts.next().unwrap().parse().unwrap();

    let hex_array = shifted_hex_to_hex_array(shifted_hex_digits);

    let integer_value = convert_hex_array_to_int(&hex_array);

    CellID {
        cell_id_prefix: integer_value,
        dataset_suffix: dataset_suffix,
    }
}

fn shifted_hex_to_hex_array(shifted_hex_digits: &str) -> Vec<i32> {
    let hex_digits: Vec<i32> = shifted_hex_digits
        .chars()
        .map(|c| (c as i32) - ('a' as i32))
        .collect();

    hex_digits
}

fn convert_hex_array_to_int(hex_array: &[i32]) -> usize {
    let hex_string: String = hex_array.iter().map(|x| format!("{:X}", x)).collect();

    let integer_value = usize::from_str_radix(&hex_string, 16).unwrap();
    integer_value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_min_x() {
        let t = Transcript::new();

        // x coordinate in 100.
        // First test should be true, since 5 < 100
        // Seoncd test should be false since 500 > 100
        assert_eq!(
            true,
            filter_transcripts(&t, 5.0, 24000.0, 0.0, 24000.0, 20.0)
        );
        assert_eq!(
            false,
            filter_transcripts(&t, 500.0, 24000.0, 0.0, 24000.0, 20.0)
        );
    }

    #[test]
    fn check_max_x() {
        let t = Transcript::new();

        // x coordinate in 100.
        // First test should be true, since 500 > 100
        // Seoncd test should be false since 5 < 100
        assert_eq!(true, filter_transcripts(&t, 0.0, 500.0, 0.0, 24000.0, 20.0));
        assert_eq!(false, filter_transcripts(&t, 0.0, 5.0, 0.0, 24000.0, 20.0));
    }

    #[test]
    fn check_min_y() {
        let t = Transcript::new();

        // y coordinate in 100.
        // First test should be true, since 5 < 100
        // Seoncd test should be false since 500 > 100
        assert_eq!(
            true,
            filter_transcripts(&t, 0.0, 24000.0, 5.0, 24000.0, 20.0)
        );
        assert_eq!(
            false,
            filter_transcripts(&t, 0.0, 24000.0, 500.0, 24000.0, 20.0)
        );
    }

    #[test]
    fn check_max_y() {
        let t = Transcript::new();

        // y coordinate in 100.
        // First test should be true, since 500 > 100
        // Seoncd test should be false since 5 < 100
        assert_eq!(true, filter_transcripts(&t, 0.0, 24000.0, 0.0, 500.0, 20.0));
        assert_eq!(false, filter_transcripts(&t, 0.0, 24000.0, 0.0, 5.0, 20.0));
    }

    #[test]
    fn check_qc() {
        let t = Transcript::new();

        // qc in 25.
        // First test should be true, since 25 > 20
        // Seoncd test should be false since 25 < 30
        assert_eq!(
            true,
            filter_transcripts(&t, 0.0, 24000.0, 0.0, 24000.0, 20.0)
        );
        assert_eq!(
            false,
            filter_transcripts(&t, 0.0, 24000.0, 0.0, 24000.0, 30.0)
        );
    }

    #[test]
    fn feature_names() {
        let mut t = Transcript::new();

        // Should be removed
        t.feature_name = "BLANK_test".to_string();
        assert_eq!(true, exclude_feature_names(&t));
        assert_eq!(
            false,
            filter_transcripts(&t, 0.0, 24000.0, 0.0, 24000.0, 20.0)
        );

        // Should be removed
        t.feature_name = "NegControlProbe_test".to_string();
        assert_eq!(true, exclude_feature_names(&t));
        assert_eq!(
            false,
            filter_transcripts(&t, 0.0, 24000.0, 0.0, 24000.0, 20.0)
        );

        // Should NOT be removed
        t.feature_name = "1_NegControlProbe_test".to_string();
        assert_eq!(false, exclude_feature_names(&t));
        assert_eq!(
            true,
            filter_transcripts(&t, 0.0, 24000.0, 0.0, 24000.0, 20.0)
        );
    }

    #[test]
    fn reverse_cell_id() {
        let cell_id = "ffkpbaba-1";

        let decoded = decode_cell_id(cell_id);

        assert_eq!(decoded.dataset_suffix, 1);
        assert_eq!(decoded.cell_id_prefix, 1437536272);
    }
}
