PACKAGE_NAME  := plotmol
CONDA_ENV_RUN := conda run --no-capture-output --name $(PACKAGE_NAME)

EXAMPLES_SKIP :=
EXAMPLES := $(filter-out $(EXAMPLES_SKIP), $(wildcard examples/*.ipynb))

.PHONY: env lint format test  test-examples

env:
	mamba create     --name $(PACKAGE_NAME)
	mamba env update --name $(PACKAGE_NAME) --file devtools/envs/base.yaml
	$(CONDA_ENV_RUN) pip install --no-build-isolation --no-deps -e .
	$(CONDA_ENV_RUN) pre-commit install || true

lint:
	$(CONDA_ENV_RUN) isort --check-only $(PACKAGE_NAME)
	$(CONDA_ENV_RUN) black --check      $(PACKAGE_NAME)
	$(CONDA_ENV_RUN) flake8             $(PACKAGE_NAME)

format:
	$(CONDA_ENV_RUN) isort  $(PACKAGE_NAME)
	$(CONDA_ENV_RUN) black  $(PACKAGE_NAME)
	$(CONDA_ENV_RUN) flake8 $(PACKAGE_NAME)

test:
	$(CONDA_ENV_RUN) pytest -v --cov=$(PACKAGE_NAME) --cov-report=xml --color=yes $(PACKAGE_NAME)/tests/

test-examples:
	$(CONDA_ENV_RUN) jupyter nbconvert --to notebook --execute $(EXAMPLES)

