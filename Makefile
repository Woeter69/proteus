.PHONY: up setup

setup:
	conda env update --file environment.yml --prune

up:
	./dev.sh
