.PHONY: up setup

setup:
	conda env update --file environment.yml --prune

up:
	./platform/dev.sh
