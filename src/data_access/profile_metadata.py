import duckdb
import pandas as pd
import numpy as np


def generate_dataframe(con, url):
    df = con.execute(
        "DESCRIBE SELECT * FROM read_parquet(?)",
        [url],
    ).fetchdf()

    df = df.drop(
        columns=["column_type", "null", "key", "default", "extra"],
        errors="ignore",
    )

    return df


def unique_counts(con, df, url, name):

    unique_counts = {}
    
    for column in df["column_name"]:
        query = f"""
            SELECT COUNT(DISTINCT "{column}")
            FROM read_parquet(?)
        """
        try:
            unique_counts[column] = con.execute(
                query, [url]
            ).fetchone()[0]
        except Exception:
            unique_counts[column] = None

    df.insert(
        1,
        f"{name}",
        df["column_name"].map(unique_counts),
    )

    return df


def concate_df(name, base_url, files):

    con = duckdb.connect()

    try:
        con.execute("INSTALL httpfs")
        con.execute("LOAD httpfs")

        dataset_url = f"{base_url.rstrip('/')}/{files[0]}"
        sample_url = f"{base_url.rstrip('/')}/{files[1]}"
        run_url = f"{base_url.rstrip('/')}/{files[2]}"

        # sample metadata
        sample_df = generate_dataframe(con, sample_url)
        sample_df = unique_counts(con, sample_df, sample_url, name)

        # dataset metadata
        dataset_df = generate_dataframe(con, dataset_url)
        dataset_df = unique_counts(con, dataset_df, dataset_url, name)

        # run metadata
        run_df = generate_dataframe(con, run_url)
        run_df = unique_counts(con, run_df, run_url, name)

        # concat to single df
        overview = pd.concat(
            [
                sample_df,
                dataset_df.iloc[[4, 6]],
                run_df.iloc[[2, 3, 4, 5, 6, 7, 8, 9]],
            ],
            ignore_index=True,
        )

        return overview

    finally:
        con.close()
        

def profile_metadata(exps):

    overviews = []

    for name, base_url, files in exps:
        overview = concate_df(name, base_url, files)
        overviews.append(overview)

    # combine df from all experiments
    combined_overview = overviews[0].iloc[:, [0]].copy()

    for overview in overviews:
        combined_overview[overview.columns[1]] = overview.iloc[:, 1]

    print(combined_overview.to_string(index=False))

    return combined_overview


'''
def profile_metadata(name, base_url, files):
    con = duckdb.connect()

    try:
        con.execute("INSTALL httpfs")
        con.execute("LOAD httpfs")

        dataset_url = f"{base_url.rstrip('/')}/{files[0]}"
        sample_url = f"{base_url.rstrip('/')}/{files[1]}"
        run_url = f"{base_url.rstrip('/')}/{files[2]}"

        # sample metadata
        sample_df = generate_dataframe(con, sample_url)
        sample_df = unique_counts(con, sample_df, sample_url, name)

        # dataset metadata
        dataset_df = generate_dataframe(con, dataset_url)
        dataset_df = unique_counts(con, dataset_df, dataset_url, name)

        # run metadata
        run_df = generate_dataframe(con, run_url)
        run_df = unique_counts(con, run_df, run_url, name)

        # concat to single df
        overview = pd.concat(
            [
                sample_df,
                dataset_df.iloc[[4, 6]],
                run_df.iloc[[2, 3, 4, 5, 6, 7, 8, 9]],
            ],
            ignore_index=True,
        )

        print(name)
        print(overview.to_string(index=False))

        return overview

    finally:
        con.close()
'''