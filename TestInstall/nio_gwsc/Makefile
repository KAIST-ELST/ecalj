include ../Makefile.define

# directroy where input/output files are saved.
TARGET=nio_gwsc
# material in the tests.
MATERIAL=nio

all: help

help:
	@echo "Do not call this Makefile alone."

# exec default case for $(TARGET).
$(TARGET).exec:
	# exec lmfa command.
	$(bindir)/lmfa $(MATERIAL) > llmfa
	# exec gwsc command.
	rm log.$(MATERIAL) QPU
	$(bindir)/gwsc 0 -np $(mpi_size) $(MATERIAL)
	@echo  >> ../summary.txt
	@-$(bindir)/dqpu QPU $(REFDIR)/QPU && echo "PASSED: $(TARGET) QPU" >> ../summary.txt || echo "FAILED: $(TARGET) QPU compared by ./bin/dqpu" >> ../summary.txt
	@-$(bindir)/diffnum  log.$(MATERIAL) $(REFDIR)/log.$(MATERIAL) 'fp evl' && echo "PASSED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt || echo "FAILED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt

# check default case for $(TARGET).
$(TARGET).check:
	# add summary.txt
	@echo  >> ../summary.txt
	@-$(bindir)/dqpu QPU $(REFDIR)/QPU && echo "PASSED: $(TARGET) QPU" >> ../summary.txt || echo "FAILED: $(TARGET) QPU compared by ./bin/dqpu" >> ../summary.txt
	@-$(bindir)/diffnum  log.$(MATERIAL) $(REFDIR)/log.$(MATERIAL) 'fp evl' && echo "PASSED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt || echo "FAILED: $(TARGET) diffnum log.$(MATERIAL) fp evl" >> ../summary.txt

# exec case1 for $(TARGET).
$(TARGET).case1.exec:
	@echo "nothing to exec for $@"

# check case1 for $(TARGET).
$(TARGET).case1.check:
	@echo "nothing to check for $@"
