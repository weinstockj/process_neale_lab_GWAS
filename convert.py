import polars as pl
import numpy as np
import subprocess, argparse, re, gzip, os, logging, sys


def log_check_error(condition: str, message: str) -> None:
    try:
        assert condition, message
    except AssertionError as e:
        logging.error(f"{e}")
        raise


def ensure_four(df: pl.DataFrame) -> pl.DataFrame:
    """ENSURE THAT CHR, BP, A1, AND A2 EXIST AS COLUMNS. MAY NEED TO RELY ON NEALE-LIKE PATTERN TO OBTAIN THESE COLUMNS. ALSO ADJUST NEALE-LIKE PATTERNS.
    df: GWAS dataframe.
    """

    # DEFINE SOME REGULAR EXPRESSIONS:
    neale_pat = "(?i)[0-9XY]{1,2}[_:\\+][0-9]+[_:\\+][AGCT]+[_:\\+][AGCT]+"
    rsid_pat = "(?i)rs[0-9]+"
    accepted_pat = "|".join([rsid_pat, neale_pat])
    snp_present = True if "SNP" in df.columns else False
    four_present = (
        True if all([x in df.columns for x in ["CHR", "BP", "A2", "A1"]]) else False
    )
    log_check_error(
        snp_present or four_present,
        f"Need SNP or all of CHR, BP, A2, and A1 to be present. GWAS columns at this point are: {df.columns}.",
    )

    # ENSURE SNP COLUMN IS PRESENT:
    if not snp_present:
        # MAKE SURE CHR COLUMN IS CLEAN THEN GENERATE SNP COLUMN:
        df = df.with_columns(pl.col("CHR").str.replace("chr", "").alias("CHR"))
        df = df.with_columns(
            pl.concat_str(
                [
                    pl.col("CHR"),
                    pl.col("BP"),
                    pl.col("A2").str.to_uppercase(),
                    pl.col("A1").str.to_uppercase(),
                ],
                separator="_",
            ).alias("SNP")
        )

    # CLEAN RSIDS (IF THEY EXIST) IN SNP COLUMN (SHOULD DO THIS BEFORE THE NEXT OPERATION TO AVOID WEIRDNESS WITH rs1283920:384923:A:G, for example):
    df = df.with_columns(
        pl.when(pl.col("SNP").str.contains(rsid_pat))
        .then(pl.col("SNP").str.extract("(?i)(rs[0-9]+)").str.to_lowercase())
        .otherwise(pl.col("SNP"))
        .alias("SNP")
    )

    # CLEAN NEALE-PATTERNS (IF THEY EXIST) IN SNP COLUMN:
    df = df.with_columns(
        pl.when(pl.col("SNP").str.contains(neale_pat))
        .then(
            pl.col("SNP")
            .str.extract("(?i)([0-9XY]{1,2}[_:\\+][0-9]+[_:\\+][AGCT]+[_:\\+][AGCT]+)")
            .str.replace_all(":|\\+", "_")
            .str.to_uppercase()
        )
        .otherwise(pl.col("SNP"))
        .alias("SNP")
    )

    # ENSURE CHR, BP, A2, AND A1 ARE PRESENT:
    if not four_present:
        # EXTRACT EACH COLUMN VALUE FROM NEALE-PATTERNS. IF RSID, PUT NULL VALUE:
        check_cols = ["CHR", "BP", "A2", "A1"]
        for i, this_col in enumerate(check_cols):
            if this_col not in df.columns:
                df = df.with_columns(
                    pl.when(pl.col("SNP").str.contains(neale_pat))
                    .then(pl.col("SNP").str.split("_").list.get(i))
                    .otherwise(None)
                    .alias(this_col)
                )
    df = df.with_columns(pl.col("A1").str.to_uppercase().alias("A1"))
    df = df.with_columns(pl.col("A2").str.to_uppercase().alias("A2"))
    df = df.with_columns(
        pl.col("CHR").str.replace("chr", "").str.to_uppercase().alias("CHR")
    )

    # SOMETIMES X AND Y CHROMOSOMES ARE REFERRED TO DIFFERENTLY:
    df = df.with_columns(
        pl.when(pl.col("CHR") == "23")
        .then(pl.lit("X"))
        .otherwise(pl.col("CHR"))
        .alias("CHR")
    )
    df = df.with_columns(
        pl.when(pl.col("CHR") == "24")
        .then(pl.lit("Y"))
        .otherwise(pl.col("CHR"))
        .alias("CHR")
    )

    # TRY TO RECOVER NON-COMPATIBLE SNP VALUES (IF THEY EXIST):
    df = df.with_columns(
        pl.when(~pl.col("SNP").str.contains(accepted_pat))
        .then(
            pl.concat_str(
                [
                    pl.col("CHR"),
                    pl.col("BP"),
                    pl.col("A2").str.to_uppercase(),
                    pl.col("A1").str.to_uppercase(),
                ],
                separator="_",
            )
        )
        .otherwise(pl.col("SNP"))
        .alias("SNP")
    )

    # COUNT HOW MANY ACCEPTIBLE PATTERNS REMAIN AND REMOVE UNACCEPTABLE LINES:
    accepted_count = df["SNP"].str.count_matches(accepted_pat).sum()
    unaccepted_count = df.shape[0] - accepted_count
    logging.info(
        f"{accepted_count}/{df.shape[0]} lines will be preserved while {unaccepted_count}/{df.shape[0]} will be dropped. SNPs must be rsIDs or be NEALE-like patterns."
    )
    df = df.filter(pl.col("SNP").str.contains(accepted_pat))

    return df


def lift_over(
    df: pl.DataFrame, chain_file: str, lift_tool: str, output: str
) -> pl.DataFrame:
    """LIFT FROM HG19 TO HG38. ASSUMES CHR, BP, A1, AND A2 EXIST AS COLUMNS. WILL RECREATE NEALE-LIKE PATTERNS AFTER LIFTOVER.
    df: GWAS dataframe.
    """

    # WRITE DATAFRAME IN BED FORMAT:
    desired_cols = df.columns
    df_bed = df.select(
        ("chr" + pl.col("CHR").cast(pl.Utf8).str.replace("(?i)^chr", "")).alias(
            "BEDCHR"
        ),  # Convert CHR to string
        (pl.col("BP").cast(pl.Int64) - 1).alias(
            "BEDBPminus1"
        ),  # Convert BP to integer before subtraction
        pl.col("BP").alias("BEDBP"),  # Keep BP as is
        pl.concat_str(
            [pl.col(col).cast(pl.Utf8) for col in df.columns], separator=","
        ).alias("BEDAllCols"),  # Ensure all columns are strings before concatenation
    )
    output_dir = os.path.dirname(os.path.abspath(output))
    output_basename = os.path.basename(output)
    to_be_lifted = os.path.join(output_dir, f"to_be_lifted_{output_basename}.bed")
    mapped = os.path.join(output_dir, f"mapped_{output_basename}.bed")
    unmapped = os.path.join(output_dir, f"unmapped_{output_basename}.bed")
    df_bed.write_csv(
        file=to_be_lifted, quote_style="never", include_header=False, separator="\t"
    )

    if not os.path.exists(to_be_lifted):
        IOError(
            f"File {to_be_lifted} does not exist. Please check the path and try again."
        )
    # LIFT OVER:
    if not os.path.exists(chain_file):
        IOError(
            f"File {chain_file} does not exist. Please check the path and try again."
        )
    subprocess.run(
        f"{lift_tool} {to_be_lifted} {chain_file} {mapped} {unmapped}", shell=True
    )

    # READ IN LIFTED BED FILE AND ADJUST TO DESIRED GWAS FORMAT:
    neale_pat = "(?i)[0-9XY]{1,2}[_:\\+][0-9]+[_:\\+][AGCT]+[_:\\+][AGCT]+"
    df = pl.read_csv(mapped, separator="\t", has_header=False)
    df.columns = ["BEDCHR", "BEDBPminus1", "BEDBP", "BEDAllCols"]
    df = df.with_columns(
        pl.col("BEDAllCols")
        .str.split_exact(by=",", n=len(desired_cols))
        .alias("split_col")
    )
    df = df.with_columns(
        [
            pl.col("split_col").struct.field(f"field_{i}").alias(name)
            for i, name in enumerate(desired_cols)
        ]
    ).drop(["split_col", "BEDAllCols", "CHR", "BP", "BEDBPminus1"])
    df = df.with_columns(
        pl.col("BEDCHR").str.replace("(?i)^chr", "").alias("CHR"),
        pl.col("BEDBP").alias("BP"),
        pl.when(pl.col("SNP").str.contains(neale_pat))
        .then(
            pl.concat_str(
                [
                    pl.col("BEDCHR").str.replace("(?i)^chr", ""),
                    pl.col("BEDBP"),
                    pl.col("A2"),
                    pl.col("A1"),
                ],
                separator="_",
            ).str.to_uppercase()
        )
        .otherwise(pl.col("SNP"))
        .alias("SNP"),
    )
    df = df.select(desired_cols)

    # REMOVE INTERMEDIATES:
    files_to_remove = [to_be_lifted, mapped, unmapped]
    for file_path in files_to_remove:
        os.remove(file_path)

    return df


def lift_rsids(df: pl.DataFrame, dbsnp_dir: str) -> pl.DataFrame:
    """GET CHR AND BP INFORMATION FOR RSIDS.

    df: GWAS dataframe (should only consist of RSIDS in SNP column).
    dbsnp_dir: Directory with parquet files that have "RSID" and "ID" columns. Example ID: chr8_60009_G_TA

    """

    dbsnp_dir = (
        re.sub("/+", "/", dbsnp_dir)
        if dbsnp_dir.endswith("/")
        else re.sub("/+", "/", dbsnp_dir) + "/"
    )

    # PREPARE VARIABLES FOR DBSNP LOOKUP:
    lookup_files = [x for x in os.listdir(dbsnp_dir) if "lookup" in x]
    accepted_chrs = [
        re.search("dbSNP_156.chr([0-9XY]+).lookup.parquet", x).group(1)
        for x in lookup_files
    ]
    mapper_df = pl.DataFrame(
        {
            "SNP": pl.Series([], dtype=pl.Utf8),
            "CHR_BP_INFO": pl.Series([], dtype=pl.Utf8),
        }
    )
    tracking_df = pl.DataFrame(
        {
            "SNP": pl.Series([], dtype=pl.Utf8),
            "CHR_BP_INFO": pl.Series([], dtype=pl.Utf8),
        }
    )
    df = df.with_columns(pl.lit(None).cast(pl.Utf8).alias("ID"))

    # SEARCH FOR VARIANTS IN LOOKUP TABLES BY CHR:
    for chromosome in accepted_chrs:
        lookup_file = f"{dbsnp_dir}dbSNP_156.chr{chromosome}.lookup.parquet"

        # FIND CHROMOSOME AND BASE PAIR INFORMATION FOR EACH RISD:
        # sub_df = df.filter(pl.col('CHR') == chromosome)
        if df.filter(pl.col("ID").is_null()).is_empty():
            continue
        lookup = pl.read_parquet(lookup_file)
        lookup = lookup.unique().unique(subset="RSID", keep="none")
        df = df.join(lookup, left_on="SNP", right_on="RSID", how="left", coalesce=True)
        df = df.with_columns(
            pl.when(pl.col("ID").is_not_null())
            .then(pl.col("ID"))
            .otherwise(pl.col("ID_right"))
            .alias("ID")
        ).drop("ID_right")
        # sub_df = sub_df.join(lookup, left_on = 'SNP', right_on = 'RSID', how = 'left', coalesce = False)
        # added_col = 'ID' if 'ID' in sub_df.columns else 'ID_right'
        # sub_df = sub_df.rename({added_col: 'CHR_BP_INFO'})

        # CHECK IF ANY OF THE PREVIOUSLY MISSED RSIDS ARE IN THIS LOOKUP FILE:
        # tracking_df = tracking_df.join(lookup, left_on = 'SNP', right_on = 'RSID', how = 'left', coalesce = False)
        # recovered_df = tracking_df.filter(~pl.col('ID').is_null()).select(['SNP', 'ID']).rename({'ID': 'CHR_BP_INFO'})
        # tracking_df = tracking_df.filter(pl.col('ID').is_null()).select(['SNP', 'ID']).rename({'ID': 'CHR_BP_INFO'})
        # nulls_df = sub_df.filter(pl.col('CHR_BP_INFO').is_null()).select(['SNP', 'CHR_BP_INFO'])
        # tracking_df = pl.concat([tracking_df, nulls_df], how = 'vertical')

        # UPDATE MAPPER DATAFRAME:
        # sub_df = sub_df.select(['SNP', 'CHR_BP_INFO'])
        # mapper_df = pl.concat([mapper_df, sub_df], how = 'vertical')

    # CONVERT TO DICTIONARY (JOINING GIVES SOME ISSUES):
    # mapper_dict = dict(zip(mapper_df['SNP'], mapper_df['CHR_BP_INFO']))
    # df = df.with_columns(pl.col('SNP').replace(mapper_dict).alias('CHR_BP_INFO'))

    # NOW UPDATE CHR AND BP INFO:
    df = df.with_columns(
        pl.when(~pl.col("ID").is_null())
        .then(pl.col("ID").str.split("_").list.get(0).str.replace("chr", ""))
        .otherwise(pl.col("CHR"))
        .alias("CHR")
    )
    df = df.with_columns(
        pl.when(~pl.col("ID").is_null())
        .then(pl.col("ID").str.split("_").list.get(1))
        .otherwise(pl.col("BP"))
        .alias("BP")
    )

    df = df.drop("ID")

    return df


def neale_wrangler(df: pl.DataFrame, dbsnp_dir: str) -> pl.DataFrame:
    """CONVERT NEALE-LIKE VARIANT NAMES (DEFINED BY neale_pat BELOW) TO RSIDS IF NECESSARY.

    df: GWAS dataframe
    dbsnp_dir: Directory with parquet files that have "RSID" and "ID" columns. Example ID: chr8_60009_G_T

    """

    dbsnp_dir = (
        re.sub("/+", "/", dbsnp_dir)
        if dbsnp_dir.endswith("/")
        else re.sub("/+", "/", dbsnp_dir) + "/"
    )
    neale_pat = "(?i)[0-9XY]{1,2}[_:\\+][0-9]+[_:\\+][AGCT]+[_:\\+][AGCT]+"

    # ACCOUNT FOR REVERSE COMPLEMENTS, SWITCHES, AND ADD 'chr' AS THE PREFIX:
    df = df.with_columns(
        pl.when(pl.col("SNP").str.contains(neale_pat))
        .then(
            "chr"
            + pl.col("SNP")
            .str.to_lowercase()
            .str.replace_all("a", "T")
            .str.replace_all("g", "C")
            .str.replace_all("c", "G")
            .str.replace_all("t", "A")
        )
        .otherwise(pl.col("SNP"))
        .alias("RC_SNP")
    )
    df = df.with_columns(
        pl.when(pl.col("SNP").str.contains(neale_pat))
        .then("chr" + pl.col("SNP"))
        .otherwise(pl.col("SNP"))
        .alias("SNP")
    )
    df = df.with_columns(
        pl.when(pl.col("SNP").str.contains(neale_pat))
        .then(
            pl.concat_str(
                [pl.col("CHR"), pl.col("BP"), pl.col("A1"), pl.col("A2")], separator="_"
            ).str.to_lowercase()
        )
        .otherwise(pl.col("SNP"))
        .alias("SNP_SWITCH")
    )
    df = df.with_columns(
        pl.when(pl.col("SNP_SWITCH").str.contains(neale_pat))
        .then(
            "chr"
            + pl.col("SNP_SWITCH")
            .str.to_lowercase()
            .str.replace_all("a", "T")
            .str.replace_all("g", "C")
            .str.replace_all("c", "G")
            .str.replace_all("t", "A")
        )
        .otherwise(pl.col("SNP_SWITCH"))
        .alias("RC_SNP_SWITCH")
    )
    df = df.with_columns(
        pl.when(pl.col("SNP_SWITCH").str.contains(neale_pat))
        .then("chr" + pl.col("SNP_SWITCH"))
        .otherwise(pl.col("SNP_SWITCH"))
        .alias("SNP_SWITCH")
    )

    # PREPARE VARIABLES FOR DBSNP LOOKUP:
    unique_chrs = df["CHR"].unique().sort()
    lookup_files = [x for x in os.listdir(dbsnp_dir) if "lookup" in x]
    accepted_chrs = [
        re.search("dbSNP_156.chr([0-9XY]+).lookup.parquet", x).group(1)
        for x in lookup_files
    ]
    unique_chrs = [x for x in unique_chrs if x in accepted_chrs]
    mapper_df = pl.DataFrame(
        {"SNP": pl.Series([], dtype=pl.Utf8), "RSID_SNP": pl.Series([], dtype=pl.Utf8)}
    )

    # SEARCH FOR VARIANTS IN LOOKUP TABLES BY CHR:
    for chromosome in unique_chrs:
        lookup_file = f"{dbsnp_dir}dbSNP_156.chr{chromosome}.lookup.parquet"

        # PERFORM JOINS FOR SNP IDs AND THE REVERSE COMPLEMENT SNP IDs:
        sub_df = df.filter(pl.col("CHR") == chromosome)
        lookup = pl.read_parquet(lookup_file)
        lookup = lookup.unique().unique(subset="ID", keep="none")
        snp_join = sub_df.join(
            lookup, left_on="SNP", right_on="ID", how="left", coalesce=False
        )
        snp_switch_join = sub_df.join(
            lookup, left_on="SNP_SWITCH", right_on="ID", how="left", coalesce=False
        )
        rc_snp_join = sub_df.join(
            lookup, left_on="RC_SNP", right_on="ID", how="left", coalesce=False
        )
        rc_snp_switch_join = sub_df.join(
            lookup, left_on="RC_SNP_SWITCH", right_on="ID", how="left", coalesce=False
        )
        added_col = "RSID" if "RSID" in snp_join.columns else "RSID_right"

        # MAKE SURE WE ACCOUNT FOR REPEAT SNP IDs IN THE LOOKUP TABLES:
        snp_join = snp_join.with_columns(
            pl.when(pl.col(added_col).n_unique().over("SNP") > 1)
            .then(None)
            .otherwise(pl.col(added_col))
            .alias("RSID")
        )
        snp_switch_join = snp_switch_join.with_columns(
            pl.when(pl.col(added_col).n_unique().over("SNP_SWITCH") > 1)
            .then(None)
            .otherwise(pl.col(added_col))
            .alias("RSID")
        )
        rc_snp_join = rc_snp_join.with_columns(
            pl.when(pl.col(added_col).n_unique().over("RC_SNP") > 1)
            .then(None)
            .otherwise(pl.col(added_col))
            .alias("RSID")
        )
        rc_snp_switch_join = rc_snp_switch_join.with_columns(
            pl.when(pl.col(added_col).n_unique().over("RC_SNP_SWITCH") > 1)
            .then(None)
            .otherwise(pl.col(added_col))
            .alias("RSID")
        )

        # COALESCE RESULTS, SUBSET TO 2 COLUMNS, AND CONCATENATE:
        sub_df = sub_df.with_columns(
            pl.when(snp_join["RSID"].is_not_null())
            .then(snp_join["RSID"])
            .when(snp_switch_join["RSID"].is_not_null())
            .then(snp_switch_join["RSID"])
            .when(rc_snp_join["RSID"].is_not_null())
            .then(rc_snp_join["RSID"])
            .when(rc_snp_switch_join["RSID"].is_not_null())
            .then(rc_snp_switch_join["RSID"])
            .otherwise(pl.col("SNP"))
            .alias("RSID_SNP")
        )

        sub_df = sub_df.select(["SNP", "RSID_SNP"])
        mapper_df = pl.concat([mapper_df, sub_df], how="vertical")

    mapper_dict = dict(zip(mapper_df["SNP"], mapper_df["RSID_SNP"]))
    df = df.with_columns(pl.col("SNP").replace(mapper_dict).alias("SNP"))
    df = df.drop(["SNP_SWITCH", "RC_SNP", "RC_SNP_SWITCH"])
    return df


def rename_col(
    new_name: str, pattern: str, df: pl.DataFrame, **kwargs
) -> pl.DataFrame | str:
    """IF DESIRED COLUMN EXISTS IN DATA FRAME, STANDARDIZE THE NAME. THROW ERROR IF NOT PRESENT AND REQUIRED.

    new_name: The name we want the column to be
    pattern: Regular expression that will match variations of the column name in the dataframe
    df: GWAS dataframe
    required: Boolan indicating if we require the column. True by default.
    return_current: Boolean indicating if we should just return the name of the given column as it currently exists. False by default.

    """

    cols = np.array(df.columns)
    required = kwargs.get("required", True)
    return_current = kwargs.get("return_current", False)
    match_list = np.array([re.search(pattern, col, re.IGNORECASE) for col in cols])
    mask_list = match_list != None

    if sum(mask_list) > 1:
        matched_names = [x.group() for x in match_list[mask_list]]
        # UNIQUE SOLUTION FOR HARMONIZED GWAS SUMMARY STATISTICS:
        if len([x for x in matched_names if x.startswith("hm_")]) == 1:
            old_name = [x for x in matched_names if x.startswith("hm_")][0]
            return old_name if return_current else df.rename({old_name: new_name})
        else:
            raise ValueError(
                f"Required pattern {pattern} does not appear exactly once in {cols}. Either pattern or columns need adjusting."
            )
    elif sum(mask_list) == 1:
        old_name = cols[mask_list].item()
        return old_name if return_current else df.rename({old_name: new_name})
    elif not required:
        return "" if return_current else df
    else:
        raise ValueError(
            f"Required pattern {pattern} does not appear exactly once in {cols}. Either pattern or columns need adjusting."
        )


def reformat_gwas(
    df: pl.DataFrame,
    dbsnp_dir: str,
    col_mapper: dict,
    initial_build: str,
    chain_file: str,
    lift_tool: str,
    output: str,
) -> pl.DataFrame:
    """REFORMAT GWAS DATAFRAME TO BE LDPRED-FUNCT COMPATIBLE.

    df: GWAS dataframe
    dbsnp_dir: Directory with parquet files that have "RSID" and "ID" columns. Example ID: chr8_60009_G_T (positions are hg38).
    col_mapper: Keys are standardized column names. Values are regular expression patterns to match those names.

    """

    # ADJUST FOR COLUMNS USED IN NEALE-LIKE NAMING CONVENTIONS FOR SNPS:
    pl.Config.set_tbl_cols(100)
    logging.info("Beginning handling of SNP, BP, A2, and A1 columns.")
    df = rename_col("SNP", col_mapper["SNP"], df, required=False)
    required = True if "SNP" not in df.columns else False
    for new_name in ["CHR", "BP", "A2", "A1"]:
        df = rename_col(new_name, col_mapper[new_name], df, required=required)
    df = ensure_four(df)
    logging.info("SNP, BP, A2, and A1 columns successfully formatted.")

    # HANDLE RSIDS DIFFERENTLY FROM NEALE-LIKE VARIANTS:
    rsid_df = df.filter(
        (pl.col("SNP").str.contains("(?i)rs[0-9]+"))
        & ~((pl.col("CHR").is_null()) | (pl.col("BP").is_null()))
    )
    rsid_null_df = df.filter(
        (pl.col("SNP").str.contains("(?i)rs[0-9]+"))
        & ((pl.col("CHR").is_null()) | (pl.col("BP").is_null()))
    )
    neale_df = df.filter(
        pl.col("SNP").str.contains(
            "(?i)[0-9XY]{1,2}[_:\\+][0-9]+[_:\\+][AGCT]+[_:\\+][AGCT]+"
        )
    )

    # BEFORE WE ALTER, STORE THE COLUMN NAMES AND TYPES:
    standard_col_names = []
    standard_col_dtypes = []
    for col_name, col_type in zip(rsid_df.columns, rsid_df.dtypes):
        standard_col_names.append(col_name)
        standard_col_dtypes.append(col_type)
    logging.info(
        f"Neale-like variant ID count: {neale_df.shape[0]}. rsID count: {rsid_df.shape[0] + rsid_null_df.shape[0]}"
    )
    logging.info(f"rsIDs with missing CHR or BP values count: {rsid_null_df.shape[0]}")

    # CHECK FOR RSIDS WITH MISSING CHR/BP VALUES. LIFTOVER TO HG38 IF NECESSARY:
    logging.info("Checking for missing CHR or BP values for RSIDs.")
    if not rsid_null_df.is_empty():
        rsid_null_df = lift_rsids(rsid_null_df, dbsnp_dir)
    logging.info("Checking for LiftOver.")
    if initial_build == "hg19":
        # HANDLE RSIDS DIFFERENTLY FROM NEALE-LIKE VARIANTS:
        if not neale_df.is_empty():
            neale_df = lift_over(neale_df, chain_file, lift_tool, output)
        # LIFTOVER IS NOT AS RELIABLE AS CONVERTING USING RSIDS MAPPED TO PARTICULAR LOCATIONS IN HG38:
        if not rsid_df.is_empty():
            rsid_df = lift_rsids(rsid_df, dbsnp_dir)

    # CONVERT NEALE-LIKE SNPS TO RSIDS IF NECESSARY:
    logging.info("Checking for Neale-like ID to RSID conversion.")
    if not neale_df.is_empty():
        neale_df = neale_wrangler(neale_df, dbsnp_dir)

    # CONCATENATE FIXED DATAFRAMES:
    for col_name, col_type in zip(standard_col_names, standard_col_dtypes):
        # print(f"Column: {col_name}, Target Type: {col_type}")
        # print(neale_df[col_name].head(5))
        neale_df = neale_df.with_columns(neale_df[col_name].cast(col_type))
        rsid_df = rsid_df.with_columns(rsid_df[col_name].cast(col_type))
        rsid_null_df = rsid_null_df.with_columns(rsid_null_df[col_name].cast(col_type))
    rsid_df = rsid_df.filter(~pl.col("SNP").is_in(neale_df["SNP"].to_list()))
    df = pl.concat([rsid_df, rsid_null_df, neale_df], how="vertical")

    # RENAME REMAINING COLUMNS:
    logging.info("Rename remaining columns.")
    for new_name, pattern in col_mapper.items():
        # DETERMINE IF COLUMN IS REQUIRED:
        if new_name in ["SNP", "CHR", "BP", "A2", "A1"]:
            continue
        elif new_name in ["N", "Z", "AF1"]:
            required = False
        else:
            required = True

        df = rename_col(new_name, col_mapper[new_name], df, required=required)

    # MAKE ADJUSTMENTS BASED ON UNIQUE CIRCUMSTANCES:
    logging.info("Creating Z column if not already present.")
    #    df = df.filter(~pl.col('BETA').str.contains('NA|^99$|^-99$'), ~pl.col('SE').str.contains('NA|^99$|^-99$'))
    df = df.filter(
        (pl.col("BETA") != 99) & (pl.col("BETA") != -99) & pl.col("BETA").is_not_null(),
        (pl.col("SE") != 99) & (pl.col("SE") != -99) & pl.col("SE").is_not_null(),
    )
    df = df.with_columns(pl.col("BETA").cast(pl.Float64), pl.col("SE").cast(pl.Float64))
    if "Z" not in df.columns:
        df = df.with_columns((pl.col("BETA") / pl.col("SE")).alias("Z"))
    logging.info(
        "Handling Neale GWAS minor_AF exception. Adding dummy AF1 column if all else fails."
    )
    if "minor_AF" in df.columns and "minor_allele" in df.columns:
        #        df = df.filter(~pl.col('minor_AF').str.contains('NA|^99$|^-99$'))
        df = df.filter(
            (pl.col("minor_AF") != 99)
            & (pl.col("minor_AF") != -99)
            & pl.col("minor_AF").is_not_null()
        )
        df = df.with_columns(pl.col("minor_AF").cast(pl.Float64))
        df = df.filter(pl.col("minor_AF") >= 0.01)
        df = df.with_columns(
            pl.when(pl.col("minor_allele") == pl.col("A1"))
            .then(pl.col("minor_AF"))
            .otherwise(1 - pl.col("minor_AF"))
            .alias("AF1")
        )
    if "AF1" not in df.columns:
        df = df.with_columns(pl.lit(0.5).alias("AF1"))

    df = df.with_columns(
        pl.when(pl.col("Z").is_infinite()).then(pl.lit(None)).otherwise(pl.col("Z"))
    )
    return df


if __name__ == "__main__":
    # DEFINE ARGUMENTS:
    parser = argparse.ArgumentParser(
        description="Reformat a GWAS summary statistics file to be pipeline compatible. Can be gzipped. Assume TSV."
    )
    parser.add_argument(
        "-s",
        "--sumstats",
        type=str,
        default="",
        help="Directory to GWAS summary statistics.",
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Name of the reformatted summary statistics."
    )
    parser.add_argument(
        "-d",
        "--dbsnp_direc",
        type=str,
        help='Directory with parquet files that have "RSID" and "ID" columns. Example ID: chr8_60009_G_T.',
    )
    parser.add_argument(
        "-ib",
        "--initial_build",
        type=str,
        default="hg38",
        help="Either hg19 or hg38. Will automatically lift over to hg38 if in hg19.",
    )
    parser.add_argument(
        "-cf",
        "--chain_file",
        type=str,
        default="/data/abattle4/seraj/tools/support_files/chain_files/hg19ToHg38.over.chain.gz",
        help="Path to hg19 to hg38 chain file.",
    )
    parser.add_argument(
        "-lt",
        "--lift_tool",
        type=str,
        default="liftOver",
        help="Path to liftOver tool.",
    )
    parser.add_argument(
        "-vf",
        "--variant_file",
        type=str,
        default="/home/jwein22/weinstocklab/projects/PRS/SNP_annotations/output/all_variants.bed",
        help="Path to variant file for final filtering.",
    )
    args = parser.parse_args()
    
    # Validate required arguments
    if not args.sumstats:
        print("Usage: python convert.py <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    # CONFIGURE LOGGING:
    log_name = re.sub(".tsv.gz", "", args.output)
    logging.basicConfig(
        filename=f"{log_name}.log",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # DEFINE COLUMN NAME MAPPINGS:
    col_mapper = {
        "SNP": "^snp(s)?$|^rs[\\-]?id(s)?$|^variant(s)?$|^id(s)?$|^marker[_\\+\\-]?name(s)?$|^name(s)?$|^hm_rsid$",
        "CHR": "^chr([_\\-\\+\\(]?(b|hg)[0-9]+\\)?)?$|^#?chrom(osome(s)?)?([_\\-\\+\\(]?(b|hg)[0-9]+\\)?)?$|^hm_chrom$",
        "BP": "^bp([_\\-\\+\\(]?(b|hg)[0-9]+\\)?)?$|^pos(ition)?([_\\-\\+\\(]?(b|hg)[0-9]+\\)?)?$|^base[_\\-\\+]?pair([_\\-\\+]?loc(ation)?)?$|^hm_pos$",
        "A2": "^a2$|^ref([_\\-]?allele)?$|^no(n|t)[_\\+\\-]?eff(ect)?[_\\+\\-]?allele$|^other[_\\+\\-]?allele$|^allele2$|^hm_other_allele$",
        "A1": "^a1$|^alt([_\\-]?allele)?$|^eff(ect)?[_\\+\\-]?allele$|^tested[_\\+\\-]?allele$|^allele1$|^hm_effect_allele$",
        "AF1": "(^(.*a1|.*alt|.*eff(ect)?|test(ed)?)[_\\+\\-]?(allele)?[_\\+\\-]?(fr(e)?q(s|uency)?)$)|(^(fr(e)?q(s|uency)?)[_\\+\\-]?(a1|alt|eff(ect)?|test(ed)?)[_\\+\\-]?(allele)?$)|^af1?$|^fr(e)?q1$|^af_alt$",
        "BETA": "^beta(s)?$|^eff(ect)?(s)?[_\\-\\+]?(size(s)?)?$|^log_odd(s)?$|^hm_beta$",
        "SE": "^se(_.*)?(beta)?$|^stand(ard)?[_\\-\\+]?err(or)?(s)?$|^log_odd(s)?_se$|^std[_\\+\\-]?err$",
        "P": "^p([_\\-\\+].*)?$|^pval(ue)?$",
        "N": "^n$|^n(_complete)?_samples$",
        "Z": "^z$|^z[_\\-\\+]?score(s)?",
    }
    accepted_builds = ["hg19", "hg38"]
    log_check_error(
        args.initial_build in accepted_builds,
        f"Initial build is {args.initial_build} but it must be hg19 or hg38.",
    )

    # IDENTIFY CHR COLUMN SO WE READ IT IN AS A STRING TO AVOID TYPE ERRORS:
    try:
        df = pl.read_csv(args.sumstats, n_rows=0, separator="\t")
        chr_col = rename_col(
            "CHR", 
            col_mapper["CHR"], 
            df, 
            return_current=True, 
            required=False
        )
        logging.info(f"CHR column has been identified as: {chr_col}")
        so = dict() if chr_col == "" else {chr_col: pl.Utf8}
        df = pl.read_csv(args.sumstats, separator="\t", schema_overrides=so)
        if "low_confidence_variant" in df.columns:
            df = df.drop(["low_confidence_variant"])
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        sys.exit(1)

    # REFORMAT GWAS AND WRITE:
    df = reformat_gwas(
        df,
        args.dbsnp_direc,
        col_mapper,
        args.initial_build,
        args.chain_file,
        args.lift_tool,
        args.output,
    )
    logging.info("Finished reformatting GWAS.")
    cols_to_select = [k for k in col_mapper.keys() if k in df.columns]
    df = df.select(cols_to_select)
    ## APK ADDED ##
    ambiguous_pairs = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

    df = (
        df.with_columns(
            pl.col("CHR").cast(pl.Int64, strict=False).alias("CHR")
         )
        .filter((pl.col("CHR") >= 1) & (pl.col("CHR") <= 22))
        .filter(
            ~pl.any_horizontal(
                pl.col("*").is_null() | pl.col("*").cast(pl.Utf8).str.contains("NaN")
            )
        )
    ).filter(
            ~(
                (pl.col("A1").cast(str) + pl.col("A2").cast(str)).is_in(
                    [f"{a}{b}" for a, b in ambiguous_pairs]
                )
            )
    )
    # Apply the awk-like transformation: adjust MAF and create final format
    df_final = df.lazy().with_columns([
            pl.concat_str([
                pl.lit("chr"),
                pl.col("CHR").cast(pl.Utf8),
                pl.lit("_"),
                pl.col("BP").cast(pl.Utf8),
                pl.lit("_"),
                pl.col("A2").cast(pl.Utf8),
                pl.lit("_"),
                pl.col("A1").cast(pl.Utf8)
            ]).alias("SNP_formatted"),
            pl.when(pl.col("AF1") < 0.5)
            .then(pl.col("AF1"))
            .otherwise(1 - pl.col("AF1"))
            .alias("MAF")
        ]).select([
            pl.col("SNP_formatted").alias("SNP"),
            pl.col("MAF"),
            pl.col("N"),
            pl.col("BETA"),
            pl.col("SE"),
            pl.col("P").alias("PVALUE")
        ])
    
    logging.info("Performing final join with variant file.")
    
    variant_df = pl.scan_csv(
                    args.variant_file, 
                    separator="\t" 
                ).select(pl.col("ID").alias("SNP"))
    
    # Filter to ensure SNPs are a subset of those in the variant file
    final_df = df_final.join(
        variant_df,
        left_on="SNP",
        right_on="SNP",
        how="inner"
    ).collect()
    
    logging.info("Writing final formatted GWAS.")
    with open(args.output, "w") as f:
        final_df.write_csv(f, separator="\t", quote_style="never")

    logging.info("Finished writing formatted GWAS.")
