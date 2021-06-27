
all:Wft Proj Lsfit Minc Depth Cwt
.PHONY:clean

SOURS1=Wft.f90 
OBJS1=Wft.o
BIN1=Wft

SOURS2=Proj.f90
OBJS2=Proj.o
BIN2=Proj

SOURS3=Lsfit.f90
OBJS3=Lsfit.o
BIN3=Lsfit


SOURS4=Minc.f90
OBJS4=Minc.o
BIN4=Minc

SOURS5=Depth.f90
OBJS5=Depth.o
BIN5=Depth

SOURS6=Cwt.f90 
OBJS6=Cwt.o
BIN6=Cwt


F90 = gfortran

${BIN1}:${OBJS1}
	${F90} -o ${BIN1} ${OBJS1}
${OBJS1}: %.o :%.F90
	${F90} -c $<

${BIN2}:${OBJS2}
	${F90} -o ${BIN2} ${OBJS2}
${OBJS2}: %.o :%.F90
	${F90} -c $<

${BIN3}:${OBJS3}
	${F90} -o ${BIN3} ${OBJS3}
${OBJS3}: %.o :%.F90
	${F90} -c $<

${BIN4}:${OBJS4}
	${F90} -o ${BIN4} ${OBJS4}
${OBJS4}: %.o :%.F90
	${F90} -c $<

${BIN5}:${OBJS5}
	${F90} -o ${BIN5} ${OBJS5}
${OBJS5}: %.o :%.F90
	${F90} -c $<

${BIN6}:${OBJS6}
	${F90} -o ${BIN6} ${OBJS6}
${OBJS6}: %.o :%.F90
	${F90} -c $<


.PHONY:clean
clean:
	-rm -f *.o 
.PHONY:uninstall
uninstall:
	-rm -f Wft Proj Lsfit Minc Depth Cwt
	
	


