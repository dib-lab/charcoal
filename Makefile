all:
	@echo "You can run the folllowing make targets: quicktest and test"

quicktest:
	py.test -k "not snakemake" tests

test:
	py.test tests
