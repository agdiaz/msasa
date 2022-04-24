.PHONY: test

# https://github.com/Nuitka/Nuitka
build:
	python -m nuitka --static-libpython=no --follow-imports --plugin-enable=numpy --plugin-enable=pyqt5 --assume-yes-for-downloads --standalone --lto src/msa.py

test:
	python ./src/msa.py \
		--input ./test/dummy.fa \
		--output ./test/dummy.afa \
		--comparer matching \
		--optimization min \
		--n-iterations 2000 \
		--output-best-plot ./test/example.png \
		--output-temp-plot ./test/example-temp.png \
		--temperature 10 \
		--engine pandas \
		--debug True