import yaml
from profile_metadata import profile_metadata

def main():
    with open("data/test_data/datasets.yaml", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    for dataset in config["datasets"]:
        profile_metadata(
            name=dataset["name"],
            base_url=dataset["base_url"],
            files=dataset["files"],
        )

if __name__ == "__main__":
    main()