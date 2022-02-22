# https://github.com/Nuitka/Nuitka
build:
	python -m nuitka --static-libpython=no --follow-imports --standalone src/msa.py

test:
	python ./src/msa.py --input ../test/dummy.fa --output ../test/example.fa --comparer matching --optimization min --n-iterations 20000 --output-best-plot ~/msasa/test/example.png --output-temp-plot ../test/example-temp.png --temperature 1 --engine numpy --debug True