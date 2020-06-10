### Makefile for GeneTack program
### $Id$

CC = gcc
LINK = g++
# http://stackoverflow.com/a/13690645/310453
D_SRC = src
OBJ = $(D_SRC)/hmm.o $(D_SRC)/getopt.o $(D_SRC)/params.o $(D_SRC)/path.o $(D_SRC)/sequence.o $(D_SRC)/fsmark.o $(D_SRC)/base_util.o
program = genetack

# $@ and $^  --  are the left and right sides of the :, respectively
# $<         --  is the first item in the dependencies list
$(D_SRC)/%.o: $(D_SRC)/%.cpp
	$(CC) -O3 -c $< -o $@

$(program): $(OBJ)
	$(LINK) $(OBJ) -o $@

# The '@' at the beginning of the line says to 'make' to not print the actual command
clean:
	@ rm -vf $(OBJ) $(program)

