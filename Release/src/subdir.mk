################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Cell.cpp \
../src/Geometry.cpp \
../src/Lattice.cpp \
../src/LocalCoords.cpp \
../src/Material.cpp \
../src/Options.cpp \
../src/Parser.cpp \
../src/Point.cpp \
../src/Quadrature.cpp \
../src/Surface.cpp \
../src/Timer.cpp \
../src/Track.cpp \
../src/TrackGenerator.cpp \
../src/Universe.cpp \
../src/log.cpp \
../src/openmoc.cpp 

OBJS += \
./src/Cell.o \
./src/Geometry.o \
./src/Lattice.o \
./src/LocalCoords.o \
./src/Material.o \
./src/Options.o \
./src/Parser.o \
./src/Point.o \
./src/Quadrature.o \
./src/Surface.o \
./src/Timer.o \
./src/Track.o \
./src/TrackGenerator.o \
./src/Universe.o \
./src/log.o \
./src/openmoc.o 

CPP_DEPS += \
./src/Cell.d \
./src/Geometry.d \
./src/Lattice.d \
./src/LocalCoords.d \
./src/Material.d \
./src/Options.d \
./src/Parser.d \
./src/Point.d \
./src/Quadrature.d \
./src/Surface.d \
./src/Timer.d \
./src/Track.d \
./src/TrackGenerator.d \
./src/Universe.d \
./src/log.d \
./src/openmoc.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include/ImageMagick -O3 -Wall -c -fmessage-length=0 -std=gnu++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


