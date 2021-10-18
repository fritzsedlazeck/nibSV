#NIMBLE_DIR?=${CURDIR}/.git/NIMBLE_DIR
#export NIMBLE_DIR
# Alternatively, use --nimbleDir:${NIMBLE_DIR} everywhere
UNAME=$(shell uname)
ifeq (${UNAME},Darwin)
	install=install_name_tool -add_rpath /opt/local/lib
else
	install=echo
endif

build:
	nim c nibsv.nim
	${install} nibsv
release:
	nim c -d:release -d:danger nibsv.nim
tests:
	@# much faster than nimble
	${MAKE} -C tests
test:
	nimble test  # uses "tests/" directory by default
install:
	nimble install -y
pretty:
	find src -name '*.nim' | xargs -L1 nimpretty --maxLineLen=1024
	find tests -name '*.nim' | xargs -L1 nimpretty --maxLineLen=1024
.PHONY: tests
