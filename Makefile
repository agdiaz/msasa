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

# MSASA
build-msasa:
	docker build --tag local/msasa .

run-msasa:
	docker run --rm -it local/msasa --help

build-sif-msasa:
	spython recipe ./Dockerfile > ./docker/msasa.def
	sudo singularity build ./docker/msasa-blosum.sif ./docker/msasa.def

# Clustal Omega
build-clustalo:
	docker build --tag local/clustalo --file ./docker/Dockerfile.clustal ./docker

run-clustalo: build-clustalo
	docker run --rm -it local/clustalo --help

build-sif-clustalo:
	spython recipe ./docker/Dockerfile.clustal > ./docker/clustalo.def
	sudo singularity build ./docker/clustalo.sif ./docker/clustalo.def

# KAlign
build-kalign:
	docker build --tag local/kalign --file ./docker/Dockerfile.kalign ./docker

run-kalign:
	docker run --rm -it local/kalign --help

build-sif-kalign:
	spython recipe ./docker/Dockerfile.kalign > ./docker/kalign.def
	sudo singularity build ./docker/kalign.sif ./docker/kalign.def

# MAFFT
build-mafft:
	docker build --tag local/mafft --file ./docker/Dockerfile.mafft ./docker

run-mafft:
	docker run --rm -it local/mafft --help

build-sif-mafft:
	spython recipe ./docker/Dockerfile.mafft > ./docker/mafft.def
	sudo singularity build ./docker/mafft.sif ./docker/mafft.def

# M.U.S.C.L.E.
build-muscle:
	docker build --tag local/muscle --file ./docker/Dockerfile.muscle ./docker

run-muscle:
	docker run --rm -it local/muscle --help

build-sif-muscle:
	spython recipe ./docker/Dockerfile.muscle > ./docker/muscle.def
	sudo singularity build ./docker/muscle.sif ./docker/muscle.def

# T-Coffee
build-tcoffee:
	docker build --tag local/tcoffee --file ./docker/Dockerfile.tcoffee ./docker

run-tcoffee:
	docker run --rm -it local/tcoffee --help

build-sif-tcoffee:
	spython recipe ./docker/Dockerfile.tcoffee > ./docker/tcoffee.def
	sudo singularity build ./docker/tcoffee.sif ./docker/tcoffee.def

# NextFlow
run-workflow:
	nextflow run ./workflows/workflow.nf -c ./workflows/nextflow.conf
