import duckdb

BASE_URL = (
    "https://browse.quantms.org/quantms/datasets/"
    "PXD011896/a516cc54b8d0"
)

FILES = [
    "openms.dataset.parquet",
    "openms.sample.parquet",
    "openms.run.parquet",
]

def truncate(value, max_length=25):
    value = str(value)
    if len(value) > max_length:
        return value[:max_length - 3] + "..."
    return value

def profile_metadata():

    con = duckdb.connect()

    try:
        con.execute("INSTALL httpfs")
        con.execute("LOAD httpfs")

        for filename in FILES:
            url = f"{BASE_URL}/{filename}"

            print(f"\n{'=' * 60}")
            print(f"File: {filename}")
            print("=" * 60)

            try:

                # Schema
                schema = con.execute(
                    "DESCRIBE SELECT * FROM read_parquet(?)",
                    [url],
                ).fetchdf()
                schema["column_type"] = schema["column_type"].apply(truncate)

                # First five rows
                preview = con.execute(
                    "SELECT * FROM read_parquet(?) LIMIT 5",
                    [url],
                ).fetchdf()

                # Combine schema and preview
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

if __name__ == "__main__":
    profile_metadata()