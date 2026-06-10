BUILD_DIR = build

TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=$(BUILD_DIR)

all: $(BUILD_DIR)/Lindhard
	./$(BUILD_DIR)/Lindhard

$(BUILD_DIR)/Lindhard: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt sources/*
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

clean:
	rm -rf $(BUILD_DIR)