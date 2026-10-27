import argparse
import json
import shutil
import uuid
import yaml
from datetime import datetime, timezone
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROMPTS_DIR = PROJECT_ROOT / "prompts"
DATASETS_FILE = PROJECT_ROOT / "data" / "test_data" / "datasets.yaml"
RUNS_DIR = PROJECT_ROOT / "runs"


def load_yaml(file_path):
    """Load a YAML configuration file."""
    with open(file_path, "r", encoding="utf-8") as file:
        return yaml.safe_load(file)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare an agentic research run."
    )

    parser.add_argument(
        "--prompt-id",
        required=True,
        help="ID of the prompt YAML file"
    )
    parser.add_argument(
        "--dataset-id",
        required=True,
        help="ID of the dataset"
    )

    args = parser.parse_args()

    prompt_id = args.prompt_id
    dataset_id = args.dataset_id

    # Load prompt configuration
    prompt_file = PROMPTS_DIR / f"{prompt_id}.yaml"

    if not prompt_file.is_file():
        parser.error(f"Prompt file not found: {prompt_file}")

    prompt_config = load_yaml(prompt_file)

    # Load dataset configuration
    if not DATASETS_FILE.is_file():
        parser.error(f"Dataset configuration not found: {DATASETS_FILE}")

    dataset_config = load_yaml(DATASETS_FILE)
    datasets = dataset_config.get("datasets", {})

    if dataset_id not in datasets:
        parser.error(f"Unknown dataset ID: {dataset_id}")

    dataset = datasets[dataset_id]
    dataset_url = dataset["url"]

    # Create a unique run directory
    timestamp = datetime.now(timezone.utc)
    run_id = (
        timestamp.strftime("%Y%m%dT%H%M%SZ")
        + "_"
        + uuid.uuid4().hex[:8]
    )

    run_dir = RUNS_DIR / prompt_id / run_id
    results_dir = run_dir / "results"
    info_dir = run_dir / "run_info"

    results_dir.mkdir(parents=True, exist_ok=False)
    info_dir.mkdir(parents=True, exist_ok=False)

    # Save copies of the input configurations
    shutil.copy2(prompt_file, info_dir / "prompt.yaml")
    shutil.copy2(DATASETS_FILE, info_dir / "datasets.yaml")

    # Save parameters and prepared inputs
    run_info = {
        "run_id": run_id,
        "timestamp_utc": timestamp.isoformat(),
        "status": "prepared_not_executed",
        "prompt_id": prompt_id,
        "dataset_id": dataset_id,
        "dataset_url": dataset_url,
        "model": None,
    }

    with open(
        info_dir / "metadata.json", "w", encoding="utf-8"
    ) as file:
        json.dump(run_info, file, indent=2, ensure_ascii=False)

    # Prepare the input for a future model call
    prepared_input = {
        "prompt": prompt_config,
        "dataset_url": dataset_url,
    }

    with open(
        info_dir / "prepared_input.json", "w", encoding="utf-8"
    ) as file:
        json.dump(prepared_input, file, indent=2, ensure_ascii=False)

    # Placeholder for future results
    (results_dir / "README.txt").write_text(
        "No AI call has been made yet.\n",
        encoding="utf-8",
    )

    print("Run prepared successfully.")
    print(f"Prompt ID: {prompt_id}")
    print(f"Dataset ID: {dataset_id}")
    print(f"Dataset URL: {dataset_url}")
    print(f"Run directory: {run_dir}")


if __name__ == "__main__":
    main()
