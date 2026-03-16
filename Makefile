.PHONY: up setup

setup:
	uv venv venv
	uv pip install -r requirements.txt

up:
	./platform/dev.sh
