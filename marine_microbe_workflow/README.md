# Marine microbe protease discovery workflow

This workflow validates the chain:

`public paired-end FASTQ -> fastp -> SPAdes -> Bakta -> BLAST+ / InterProScan evidence -> OpenFold3`

Validation sample: *Alteromonas macleodii* strain SSi026, isolated from *Sargassum ilicifolium* at St John's Island, Singapore. The sample is public under ENA `SAMN47740168`; reads are `SRR32944547`; the NCBI reference assembly is `GCF_049379645.1`.

The reference annotations are used only as a held-out audit target. Candidate discovery is based on the assembled/annotated outputs and a fixed protein homology database; reference labels are revealed in the audit step.

The production deployment target is Linux/HPC or cloud. The macOS run is an execution and integration check; the Bakta light database and InterPro Matches API are explicitly marked as validation modes, not substitutes for full production databases.
