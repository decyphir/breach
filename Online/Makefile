MEX         =  /Applications/MATLAB_R2015aSP1.app/bin/mex
CXXFLAGS = #-DDEBUG__  #-DDEBUG___  #-std=c++11 -O3 #-DYYDEBUG=1  

SRC_DIR     = src/
BIN_DIR     = bin/
OBJS_DIR     = obj/
INCLUDE_DIR = include/
INCLUDE     = $(addprefix -I, $(INCLUDE_DIR))

SL           = onlineMonitorWrapper.cpp
SL_SRC       = $(SRC_DIR)$(SL)  
SL_MEX       =  $(addprefix $(BIN_DIR), $(patsubst %.cpp,%.mexmaci64,$(SL)))

CL           = stl_eval_mex.cpp
CL_SRC       = $(SRC_DIR)$(CL)  
CL_MEX       = $(addprefix $(BIN_DIR), $(patsubst %.cpp,%.mexmaci64,$(CL)))

CLT           = test_online_mex.cpp
CLT_SRC       = $(SRC_DIR)$(CLT)  
CLT_MEX       = $(addprefix $(BIN_DIR), $(patsubst %.cpp,%.mexmaci64,$(CLT)))

SRCS_FILES    = stl_parser.cpp \
          stl_scanner.cpp \
          stl_driver.cpp \
          tools.cpp \
          stl_atom.cpp \
          transducer.cpp \
          interval_transducer.cpp \
          interval.cpp \
          update_transducer.cpp \
          robustness.cpp \
          signal.cpp \
          signal_expr.cpp 

HEADRS_FILES = FlexLexer.h \
		interval.h \
		location.hh \
		position.hh \
		robustness.h \
		signal.h \
		signal_expr \
		stack.hh \
		stl_driver.h \
		stl_parser.h \
		stl_scanner.h \
		tools.h \
		transducer.h 


SRCS   = $(addprefix $(SRC_DIR), $(SRCS_FILES))
HEADRS   = $(addprefix $(INCLUDE_DIR), $(HEADRS_FILES))

OBJS_FILES  = $(patsubst %.cpp,%.o,$(SRCS_FILES))
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_FILES))

all: $(SL_MEX) $(CL_MEX) $(CLT_MEX)

obj/%.o : src/%.cpp 
	$(MEX) -c $(CXXFLAGS) $(INCLUDE) $< -outdir $(OBJS_DIR)

$(SL_MEX) : $(SL_SRC) $(OBJS)
	$(MEX) $(SL_SRC) $(CXXFLAGS) $(INCLUDE) $(OBJS) -outdir $(BIN_DIR) 

$(CL_MEX) : $(CL_SRC) $(OBJS)
	$(MEX) $(CL_SRC) $(CXXFLAGS) $(INCLUDE) $(OBJS) -outdir $(BIN_DIR) 

$(CLT_MEX) : $(CLT_SRC) $(OBJS)
	$(MEX) $(CLT_SRC) $(CXXFLAGS) $(INCLUDE) $(OBJS) -outdir $(BIN_DIR) 

clean:  
	rm -f $(OBJS) $(SL_MEX) $(CL_MEX) $(CLT_MEX)

echo:
	@echo $(HEADRS)
	@echo $(SRCS)
	@echo $(INCLUDE) 
	@echo $(OBJS)
	@echo $(SL_MEX)

