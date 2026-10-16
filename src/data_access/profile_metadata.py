import duckdb

def truncate(value, max_length=25):
    value = str(value)
    if len(value) > max_length:
        return value[:max_length - 3] + "..."
    return value

def profile_metadata(name, base_url, files):

    con = duckdb.connect()

    try:
        con.execute("INSTALL httpfs")
        con.execute("LOAD httpfs")

        for filename in files:
            url = f"{base_url.rstrip('/')}/{filename}"

            print(f"\n{'=' * 60}")
            print(f"File: {filename}")
            print("=" * 60)

            try:

                # schema
                schema = con.execute(
                    "DESCRIBE SELECT * FROM read_parquet(?)",
                    [url],
                ).fetchdf()
                schema = schema.drop(
                    columns=["column_type", "null", "key", "default", "extra"],
                    errors="ignore",
                )

                # first five rows as preview
                preview = con.execute(
                    "SELECT * FROM read_parquet(?) LIMIT 5",
                    [url],
                ).fetchdf()

                # unique values per column
                unique_counts = {}

                for column in schema["column_name"]:
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

                schema.insert(
                    1,
                    "unique_values",
                    schema["column_name"].map(unique_counts),
                )

                # combine schema and preview
                for i in range(len(preview)):
                    schema[f"row_{i + 1}"] = [
                        truncate(value) for value in preview.iloc[i].tolist()
                    ]

                print("\nSchema and preview:")
                print(schema.to_string(index=False))

            except Exception as e:
                print(f"Could not read file: {e}")

    finally:
        con.close()