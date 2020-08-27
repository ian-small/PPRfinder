PPRfinder.jar: compile
	cd src && jar cvmf ../META-INF/MANIFEST.MF ../$@  pprfinder/*.class resources

compile:
	cd src && javac -Xlint:unchecked pprfinder/*.java 

