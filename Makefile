.PHONY: up setup test

setup:
	uv venv venv
	VIRTUAL_ENV=venv uv pip install -r requirements.txt

test:
	PYTHONPATH=. ./venv/bin/pytest tests/

docs:
	./venv/bin/mkdocs serve

up:
	./platform/dev.sh
