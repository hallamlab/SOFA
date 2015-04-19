LAST = LAST
FragGeneScanPlus = FragGeneScanPlus
FLASH = FLASH

all:
	$(MAKE) -C $(LAST)	
	$(MAKE) -C $(FragGeneScanPlus)
	$(MAKE) -C $(FLASH)
