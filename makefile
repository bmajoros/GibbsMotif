CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl
#---------------------------------------------------------
$(OBJ)/sim.o:\
		sim.C
	$(CC) $(CFLAGS) -o $(OBJ)/sim.o -c \
		sim.C
#---------------------------------------------------------
sim: \
		$(OBJ)/Motif.o \
		$(OBJ)/sim.o
	$(CC) $(LDFLAGS) -o sim \
		$(OBJ)/Motif.o \
		$(OBJ)/sim.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/Motif.o:\
		Motif.C\
		Motif.H
	$(CC) $(CFLAGS) -o $(OBJ)/Motif.o -c \
		Motif.C
#---------------------------------------------------------
$(OBJ)/gibbs.o:\
		gibbs.C
	$(CC) $(CFLAGS) -o $(OBJ)/gibbs.o -c \
		gibbs.C
#---------------------------------------------------------
gibbs: \
		$(OBJ)/Motif.o \
		$(OBJ)/gibbs.o
	$(CC) $(LDFLAGS) -o gibbs \
		$(OBJ)/Motif.o \
		$(OBJ)/gibbs.o \
		$(LIBS)
#---------------------------------------------
