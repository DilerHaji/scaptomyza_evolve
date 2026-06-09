import polars as pl
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    header_str = "chrom,start,end," \
                 "slope_divergence,z_divergence,p_divergence," \
                 "slope_pbs,z_pbs,p_pbs," \
                 "slope_pbs_ref,z_pbs_ref,p_pbs_ref," \
                 "slope_pbsn1,z_pbsn1,p_pbsn1," \
                 "slope_pbsn1_ref,z_pbsn1_ref,p_pbsn1_ref," \
                 "slope_pbe,z_pbe,p_pbe," \
                 "slope_pbe_ref,z_pbe_ref,p_pbe_ref\n"


    schema_map = {
        "chrom": pl.String,
        "start": pl.Int64,
        "end": pl.Int64,
        "slope_divergence": pl.Float64, "z_divergence": pl.Float64, "p_divergence": pl.Float64,
        "slope_pbs": pl.Float64, "z_pbs": pl.Float64, "p_pbs": pl.Float64,
        "slope_pbs_ref": pl.Float64, "z_pbs_ref": pl.Float64, "p_pbs_ref": pl.Float64,
        "slope_pbsn1": pl.Float64, "z_pbsn1": pl.Float64, "p_pbsn1": pl.Float64,
        "slope_pbsn1_ref": pl.Float64, "z_pbsn1_ref": pl.Float64, "p_pbsn1_ref": pl.Float64,
        "slope_pbe": pl.Float64, "z_pbe": pl.Float64, "p_pbe": pl.Float64,
        "slope_pbe_ref": pl.Float64, "z_pbe_ref": pl.Float64, "p_pbe_ref": pl.Float64
    }

    try:
        lf = pl.concat([
            pl.scan_csv(f, schema_overrides=schema_map) 
            for f in args.inputs
        ])        
        lf.collect().write_csv(args.output)
        
    except Exception as e:
        with open(args.output, "w") as f:
            f.write(header_str)

if __name__ == "__main__":
    main()