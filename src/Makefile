$(shell mkdir -p ../obj)
$(shell mkdir -p ../bin)

export SHELL=/bin/bash
export ECHO=/bin/echo
export SRC_DIR := $(shell pwd)
export OBJ_DIR := $(subst src,obj,$(SRC_DIR))
export BIN_DIR := $(subst src,bin,$(SRC_DIR))
export COMMON_DIR := $(SRC_DIR)/SR_Common
export INCLUDES := -I$(COMMON_DIR)

export CC = gcc
export CFLAGS = -Wall -g -std=gnu99

VPATH := $(SRC_DIR)/SR_Build:$(SRC_DIR)/SR_Common

BUILD_OBJ := SR_Build_Main.o SR_Build_GetOpt.o SR_OutHashTable.o SR_Error.o SR_Reference.o md5.o
SR_BUILD_OBJ := $(addprefix $(OBJ_DIR)/,$(BUILD_OBJ))

DEP = $(BUILD_OBJ:.o=.d)
SR_BUILD_DEP = $(addprefix $(OBJ_DIR)/,$(DEP))

SUBDIR = SR_Build SR_Common

all: dep SR_Build

dep:
	@for dir in $(SUBDIR); do \
	    $(MAKE) --no-print-directory -C $$dir; \
	    echo ""; \
	done

SR_Build: $(BUILD_OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -g -o $(BIN_DIR)/SR_Build $(SR_BUILD_OBJ)
	@$(ECHO) -e "\n"


-include $(SR_BUILD_DEP)


.PHONY: SR_Build
.PHONY: all
.PHONY: dep
.PHONY: clean

clean:
	-rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

