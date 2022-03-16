SRC := src
OBJ := obj
OUT := bin

CFLAGS := -D_VERBOSITY=4

$(OBJ)/%.o: $(SRC)/%.c
	gcc -c $^ -o $@ $(CFLAGS)
	
mrhs: $(OBJ)/mrhs.bm.o $(OBJ)/mrhs.bv.o $(OBJ)/mrhs.o $(OBJ)/mrhs.hillc.o $(OBJ)/mrhs.rz.o $(OBJ)/mrhs.tester.o $(OBJ)/mrhs.1.7.o
	gcc $^ -o $(OUT)/mrhs -lm

clean:
	rm ./$(OUT)/mrhs 
	rm ./$(OBJ)/* 
 