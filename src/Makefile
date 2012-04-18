$(shell mkdir -p ../obj)
$(shell mkdir -p ../bin)

export SHELL=/bin/bash
export ECHO=/bin/echo
export SRC_DIR := $(shell pwd)
export OBJ_DIR := $(subst src,obj,$(SRC_DIR))
export BIN_DIR := $(subst src,bin,$(SRC_DIR))
export SAMTOOLS_DIR := $(SRC_DIR)/samtools
export COMMON_DIR := $(SRC_DIR)/TGM_Common
export INCLUDES := -I$(COMMON_DIR)
export LIB := -L$(OBJ_DIR) -lbam -lz -lm

export CC = gcc
export CFLAGS = -Wall -g -std=gnu99

VPATH := $(SRC_DIR)/TGM_Scan:$(SRC_DIR)/TGM_Common

SCAN_OBJ := TGM_ReadPairScanMain.o TGM_ReadPairScan.o TGM_ReadPairScanGetOpt.o TGM_Error.o TGM_BamInStream.o TGM_LibInfo.o TGM_Utilities.o TGM_FragLenHist.o \
	    TGM_BamMemPool.o TGM_BamPairAux.o TGM_BamHeader.o TGM_GetOpt.o

TGM_SCAN_OBJ := $(addprefix $(OBJ_DIR)/,$(SCAN_OBJ))

DEP = $(SCAN_OBJ:.o=.d)
TGM_SCAN_DEP = $(addprefix $(OBJ_DIR)/,$(DEP))

SUBDIR = TGM_Scan TGM_Common

all: dep samtools TGM_Scan

dep:
	@for dir in $(SUBDIR); do \
	    $(MAKE) --no-print-directory -C $$dir; \
	    echo ""; \
	done

samtools:
	make lib -C $(SAMTOOLS_DIR)
	mv $(SAMTOOLS_DIR)/libbam.a $(OBJ_DIR)/

TGM_Scan: $(SCAN_OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -g -o $(BIN_DIR)/TGM_Scan $(TGM_SCAN_OBJ) $(LIB) 
	@$(ECHO) -e "\n"


-include $(TGM_SCAN_DEP)


.PHONY: TGM_Scan
.PHONY: all
.PHONY: dep
.PHONY: clean
.PHONY: samtools

clean:
	-rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
	make clean -C $(SAMTOOLS_DIR)

