#NIMBLE_DIR?=${CURDIR}/nimbleDir
#export NIMBLE_DIR
# Alternatively, use --nimbleDir:${NIMBLE_DIR} everywhere
UNAME=$(shell uname)
ifeq (${UNAME},Darwin)
	install=install_name_tool -add_rpath /opt/local/lib
else
	install=echo
endif

build:
	nim c src/nibsv.nim
	${install} src/nibsv
release:
	nim c -d:release -d:danger src/nibsv.nim
all:
	${MAKE} install
quick:
	nim c -r tests/t_kmers.nim
	nim c -r tests/t_util.nim
help:
	nimble -h
	nimble tasks
tests:
	@# much faster than nimble
	${MAKE} -C tests
test:
	nimble test  # uses "tests/" directory by default
integ-test:
	@echo 'integ-test TBD'
install:
	nimble install -y
pretty:
	find src -name '*.nim' | xargs -L1 nimpretty --maxLineLen=1024
	find tests -name '*.nim' | xargs -L1 nimpretty --maxLineLen=1024
setup: #vendor/threadpools vendor/STRling
	nimble install --verbose -y hts kmer bitvector cligen msgpack4nim
	#cd vendor/threadpools; nimble install --verbose -y
	#cd vendor/STRling; nimble install --verbose -y
vendor/threadpools vendor/STRling:
	git submodule update --init
rsync: # not used for now
	mkdir -p ${NIMBLE_DIR}/pkgs/
	rsync -av vendor/STRling/ ${NIMBLE_DIR}/pkgs/strling-0.3.0/
	rsync -av vendor/threadpools/ ${NIMBLE_DIR}/pkgs/threadpools-0.1.0/

.PHONY: test tests
