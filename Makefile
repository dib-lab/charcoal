all:
	@echo "You can run the following make targets: quicktest and test"

quicktest:
	py.test -k "not snakemake" tests

test:
	py.test tests
