-include .build/fileIO.d
-include .build/ncoord.d
-include .build/windowInfo.d
-include .build/util.d
-include .build/testBoosting.d
-include .build/classifier.d

.build/testBoosting.l: \
	.build/fileIO.d \
	.build/ncoord.d \
	.build/windowInfo.d \
	.build/util.d \
	.build/testBoosting.d \
	.build/classifier.d \

bin/testBoosting: .build/testBoosting.l \
	.build/fileIO.o \
	.build/ncoord.o \
	.build/windowInfo.o \
	.build/util.o \
	.build/testBoosting.o \
	.build/classifier.o \

