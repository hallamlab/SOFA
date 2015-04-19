LAST = LAST
FragGeneScanPlus = FragGeneScanPlus
FLASH = FLASH

all:
	$(MAKE) -C $(LAST)	
	$(MAKE) -C $(FragGeneScanPlus)
	$(MAKE) -C $(FLASH)

clean:
	$(MAKE) -C $(LAST) clean
	$(MAKE) -C $(FragGeneScanPlus) clean
	$(MAKE) -C $(FLASH) clean

