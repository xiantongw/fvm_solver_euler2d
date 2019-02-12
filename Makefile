LIBS_DIR = ./libs
SRC_DIR = ./src
CC = icc

build_libs:
	@if([ ! -d "${LIBS_DIR}" ]); then mkdir "${LIBS_DIR}"; fi;
	$(CC) ${SRC_DIR}/flux.c -fPIC -shared -o ${LIBS_DIR}/flux.so
	$(CC) ${SRC_DIR}/euler_steady_solver.c -fPIC -shared -o ${LIBS_DIR}/euler_steady_solver.so
	
clean:
	rm -f "${LIBS_DIR}/flux.so"
	rm -f "${LIBS_DIR}/euler_steady_solver.so"
	rm -rf "./wrapper/__pycache__"
	rm -rf "./mesh/__pycache__"

clean~:
	rm -rf *~
