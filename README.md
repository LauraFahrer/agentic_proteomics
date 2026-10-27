# Agentic Proteomics
This exploratory research project investigates the capabilities and limitiations of LLMs, such as Claude, in hypothesis generation and agentic research using large publicly available proteomics datasets.

## Goals
- Explore different research questions, data sources and analytical workflows
- Investigate capabilities as well as limitations of agents
- Evaluate generated hypotheses with regard to biological plausibility, novelty and existing data support (see Evaluation Criterias)
- Integrate literature research and other omics technologies
- Document dataset and additional sources of evidence like scientific publications or databases
- Test promising hypotheses on new data or experimentally

## Tools
Tools and technologies are currently under exploration and will be documented as the project develops.

## Datasets
- Publicly available proteomics datasets of human cancer cell lines on quantms

## Evaluation Criterias
- Comprehensibility: essential statements and results are understandable and supported
- Novelty: generated hypotheses are new and go beyond known relations
- Data-driven: observations and statements are backed by the data
- Biological plausibility: generated hypotheses are scientificly coherent
- Testability: statements and hypotheses can be validated experimentally or on new data

## Running a LLM query
The src/run.py script prepares a LLM run based on a specified prompt and dataset. It loads the corresponding YAML configurations and creates a dedicated directory for documenting the run.

Current status: The script is a preliminary implementation. It prepares the run configuration but does not yet call an LLM or execute research tasks.

Usage: Run the script from the project root directory:
    python src/run.py --promt-id <PROMPT_ID> --dataset-id <DATASET_ID>

Parameters:

- `--prompt-id`: Identifier of the prompt configuration. The corresponding YAML file must exist in prompts/ and follow the naming convention <PROMPT_ID>.yaml.

- `--dataset-id`: Identifier of the dataset defined in data/test_data/datasets.yaml. The dataset configuration must include a URL.

Output: Each execution creates a unique run directory under runs/<PROMPT_ID>/.

    runs/
    └── <PROMPT_ID>/
        └── <RUN_ID>/
            ├── results/
            └── run_info/
                └── run_info.json

   - results/: Contains the research results and generated outputs
   - run_info/: Contains information about the run in a .json file including run ID timestamp, prompt ID, dataset ID, dataset URL and model. For each executed run a unique run ID is generated to prevent overwriting previous runs.