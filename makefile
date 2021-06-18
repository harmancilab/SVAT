all: SVAT

CC = g++
comp_flags = -c -O3 -Wall
lib_flags = -lz
exec_name = bin/SVAT
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

%.o: %.c
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/main.o \
${LIB_DIR}/svat_indel_utils.o \
${LIB_DIR}/svat_snv_utils.o \
${LIB_DIR}/svat_aggregation_utils.o \
${LIB_DIR}/svat_utils.o \
${LIB_DIR}/svat_secure_aggregation_utils.o \
${LIB_DIR}/svat_secure_snv_utils.o \
${LIB_DIR}/svat_ansi_string.o \
${LIB_DIR}/svat_annot_region_tools.o \
${LIB_DIR}/svat_gff_utils.o \
${LIB_DIR}/svat_ansi_cli.o \
${LIB_DIR}/svat_config.o \
${LIB_DIR}/svat_genome_sequence_tools.o \
${LIB_DIR}/svat_variation_tools.o \
${LIB_DIR}/svat_rng.o \
${LIB_DIR}/svat_seed_manager.o \
${LIB_DIR}/svat_nomenclature.o \
${LIB_DIR}/svat_nucleotide.o \
${LIB_DIR}/file_utils.o

SVAT: ${objs}
	${CC} -O3 ${lib_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 

