UNAME=$(shell uname)
ifeq (${UNAME},Darwin)
	install=install_name_tool -add_rpath /opt/local/lib
else
	install=echo
endif

tests:
	nim c all.nim
	${install} all
	./all
clean:
	rm -f all
