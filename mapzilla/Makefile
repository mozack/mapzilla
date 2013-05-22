# Make file for Mapzilla

SRCDIR=src/main/c

all: clean java

java:
	mvn clean package

zilla:
	g++ -g -pthread -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -fPIC $(SRCDIR)/mapzilla.c -o mapz

jr: jrclean jrsrc

jrclean:
	rm jr

jrsrc:
	g++ -g -pthread $(SRCDIR)/jiffyread.c -o jr

fw: fwclean fwsrc

fwclean:
	rm fw	

fwsrc:
	g++ -g -pthread $(SRCDIR)/fastwrite2.c -o fw

clean:
	mvn clean
	rm -rf target
