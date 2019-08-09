# Cattle Assembly Scripts
---

Accessory scripts used in the assembly and validation of the Cattle Reference Genome.

Some analysis depends on [GetMaskBedFasta](https://github.com/njdbickhart/GetMaskBedFasta) which is a java program that was written to rapidly identify "N" masked regions of an assembly fasta. You can run this program as a jar file as follows:

```bash
java -Xmx10G -jar GetMaskBedFasta/store/GetMaskBedFasta.jar -h
```

Without options, the program will try to create a GUI so that you can select files from a file browser. If you want to use the GUI instead of entering command line commands, then remove the "-h" after the jar file invocation.
