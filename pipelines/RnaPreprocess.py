import os

import bioluigi.tasks
from bioluigi.utils import get_ext, parseHISATLog, CheckTargetNonEmpty
from bioluigi.slurm import SlurmExecutableTask
from bioluigi.decorators import requires, inherits

import luigi
from luigi import LocalTarget, Parameter


luigi.auto_namespace(scope=__name__)


timecourse_libs = ['ERR1224553', 'ERR1224554', 'ERR1224555', 'ERR1224556', 'ERR1224557',
                   'ERR1224558', 'ERR1224559', 'ERR1224560', 'ERR1224561', 'ERR1224562',
                   'ERR1224563', 'ERR1224564', 'ERR1224565', 'ERR1224566', 'ERR1224567',
                   'ERR1224568', 'ERR1224569', 'ERR1224570', 'ERR1224571', 'ERR1224572',
                   'ERR1224573', 'ERR1224574', 'ERR1224575', 'ERR1224576', 'ERR1224577',
                   'ERR1224578', 'ERR1224579', 'ERR1224580', 'ERR1224581', 'ERR1224582',
                   'ERR1224583', 'ERR1224584', 'ERR1224585', 'ERR1224586', 'ERR1224587',
                   'ERR1224588', 'ERR1224589', 'ERR1224590', 'ERR1224591', 'ERR1224592',
                   'ERR1224593', 'ERR1224594']


PST104E137_rnaseq_libs = ['SRR6043965', 'SRR6043966', 'SRR6043967', 'SRR6043968', 'SRR6043969',
                          'SRR6043970', 'SRR6043971', 'SRR6043972', 'SRR6043973', 'SRR6043974',
                          'SRR6043975', 'SRR6043976', 'SRR6043977', 'SRR6043978', 'SRR6043979']

lib_path = '/nbi/Research-Groups/JIC/Diane-Saunders/YR_2018_genome_annotation/timecourse_rnaseq_data'
sra_path = '/nbi/Research-Groups/JIC/Diane-Saunders/YR_2018_genome_annotation/104E137_rnaseq_data'
reads_dir = '/jic/scratch/groups/Diane-Saunders/buntingd/reads'

# --------- Collecting input reads -----------#
'''Pulls in the timecourse reads a fastq from reads_dir and the PST104 paper
   RNASeq from the SRA dumps in sra_path, then runs Trimmomatic for adaptor clipping'''


class GetSRA(luigi.ExternalTask):
    library = Parameter()

    def output(self):
        return LocalTarget(os.path.join(sra_path, self.library, self.library + '.sra'))


@requires(GetSRA)
class FastqDump(bioluigi.tasks.sra.FastqDump):
    scratch_dir = Parameter()

    def output(self):
        return [LocalTarget(os.path.join(reads_dir, self.library, self.library + "_1.fastq.gz")),
                LocalTarget(os.path.join(reads_dir, self.library, self.library + "_2.fastq.gz"))]


class GetTCFastq(luigi.ExternalTask):
    library = Parameter()

    def output(self):
        return [LocalTarget(os.path.join(lib_path, self.library + '_1.fastq')),
                LocalTarget(os.path.join(lib_path, self.library + '_2.fastq'))]


@inherits(FastqDump, GetTCFastq)
class GetFastq(luigi.WrapperTask):
    base_dir = luigi.Parameter()

    def requires(self):
        if self.library in timecourse_libs:
            return self.clone(GetTCFastq)
        elif self.library in PST104E137_rnaseq_libs:
            return self.clone(FastqDump)

    def output(self):
        return self.input()


@requires(GetFastq)
class Trimmomatic(bioluigi.tasks.trimmomatic.Trimmomatic):
    def output(self):
        return [LocalTarget(os.path.join(reads_dir, self.library, self.library + "_filtered_R1.fastq.gz")),
                LocalTarget(os.path.join(reads_dir, self.library, self.library + "_filtered_R2.fastq.gz")),
                LocalTarget(os.path.join(self.base_dir, self.library, self.library + "_trimmomatic.txt"))]


# --------- Alignment to reference -----------#
'''Gets the reference and makes a HISAT index then aligns reads and sorts the bam.
   The adaptor task is necessary as HISAT produces two outputs, the bam and the alignment
   log and for downstream tasks we just need the bam.
'''


class ReferenceFasta(luigi.ExternalTask):
    reference = luigi.Parameter()

    def output(self):
        return LocalTarget(self.reference)


@requires(ReferenceFasta)
class HISATIndexGenome(bioluigi.tasks.hisat.HISATIndexGenome):
    pass


@requires(reads=Trimmomatic, genome=HISATIndexGenome)
class HISAT(bioluigi.tasks.hisat.HISAT):
    base_dir = luigi.Parameter()

    def output(self):
        return {'bam': LocalTarget(os.path.join(self.scratch_dir, get_ext(os.path.basename(self.reference))[0],
                                                self.library, self.library + '.bam')),
                'log': LocalTarget(os.path.join(self.base_dir, self.library, self.library + '.hisat.log'))}


@requires(HISAT)
class Adaptor(luigi.WrapperTask):
    def output(self):
        return self.input()['bam']


@requires(Adaptor)
class SamtoolsSort(bioluigi.tasks.samtools.SamtoolsSort):
    pass


# --------- Stringtie Transcript Reconstruction -----------#
'''Runs Stringtie on each individual bam file, then merges the output GTFs into
   on file which is then converted to GFF3.
   The merge task filters out libraries with less that 10% alignment rate from
   downstream analysis.
'''


@requires(SamtoolsSort)
class StringTie(bioluigi.tasks.stringtie.StringTie):
    pass


@inherits(StringTie, HISAT)
class StringTieMerge(bioluigi.tasks.stringtie.StringTieMerge):
    library = None

    def requires(self):
        return [self.clone(StringTie, library=lib) for lib in timecourse_libs + PST104E137_rnaseq_libs]

    def input(self):
        logs = [parseHISATLog(self.clone(HISAT, library=lib).output()['log'].path, lib)
                for lib in timecourse_libs + PST104E137_rnaseq_libs]
        filtered = [x['Library'] for x in logs if x['Overall alignment rate'] > 0.1]
        print([(x['Library'], x['Overall alignment rate']) for x in logs if x['Overall alignment rate'] > 10])
        return [self.clone(StringTie, library=lib).output() for lib in filtered]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'stringtie.gtf'))


@requires(StringTieMerge)
class StringTieGFF3(bioluigi.tasks.misc.GTFtoGFF3):
    pass


# --------- Cufflinks Transcript Reconstruction -----------#
'''Runs Cufflinks on each individual bam file, then merges the output GTFs into
   on file which is then converted to GFF3.
   The merge task filters out libraries with less that 10% alignment rate from
   downstream analysis.
'''


@requires(SamtoolsSort)
class Cufflinks(bioluigi.tasks.cufflinks.Cufflinks):
    pass


@inherits(Cufflinks, HISAT)
class CuffMerge(bioluigi.tasks.cufflinks.CuffMerge):
    library = None

    def requires(self):
        return [self.clone(Cufflinks, library=lib) for lib in timecourse_libs + PST104E137_rnaseq_libs]

    def input(self):
        logs = [parseHISATLog(self.clone(HISAT, library=lib).output()['log'].path, lib)
                for lib in timecourse_libs + PST104E137_rnaseq_libs]
        filtered = [x['Library'] for x in logs if x['Overall alignment rate'] > 0.1]
        print([(x['Library'], x['Overall alignment rate']) for x in logs if x['Overall alignment rate'] > 10])
        return [self.clone(Cufflinks, library=lib).output() for lib in filtered]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'cufflinks.gtf'))


@requires(CuffMerge)
class AddTranscripts(bioluigi.tasks.cufflinks.AddTranscripts):
    pass


@requires(AddTranscripts)
class CufflinksGFF3(bioluigi.tasks.misc.GTFtoGFF3):
    pass


# --------- Trinity Genome Guided Assembly -----------#
'''Merges all the individual bams into one then runs Trinity in
   genome guided mode.
'''


@inherits(SamtoolsSort, HISAT)
class MergeBam(bioluigi.tasks.trinity.MergeBam):
    library = None

    def requires(self):
        return [self.clone(SamtoolsSort, library=lib) for lib in timecourse_libs + PST104E137_rnaseq_libs]

    def input(self):
        logs = [parseHISATLog(self.clone(HISAT, library=lib).output()['log'].path, lib)
                for lib in timecourse_libs + PST104E137_rnaseq_libs]
        filtered = [x['Library'] for x in logs if x['Overall alignment rate'] > 0.1]
        print([(x['Library'], x['Overall alignment rate']) for x in logs if x['Overall alignment rate'] > 10])
        return [self.clone(SamtoolsSort, library=lib).output() for lib in filtered]

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, get_ext(os.path.basename(self.reference))[0], 'merged.bam'))

    # def run(self):
    #     super().run()
    #     yield self.clone((requires(self.__cls__)(bioluigi.tasks.SamtoolsIndex)))


@requires(MergeBam)
class TrinityGG(bioluigi.tasks.trinity.TrinityGG):
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'trinityGG.fasta'))


# --------- Trinity De Novo Assembly -----------#
'''Converts the merged bam back into FASTQ format by first name-sorting
   (so pair reads are next to each other) then running BAMtoFASTQ.

    Runs Trinity in de novo mode on the resulting reads.
    This is necessary to filter out host reads that would otherwise
    contaminate the assembly.
'''


@requires(Adaptor)
class SamtoolsSortName(bioluigi.tasks.samtools.SamtoolsSortName):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "nbi-medium  --tmp 20G"


@requires(SamtoolsSortName)
class BAMtoFASTQ(CheckTargetNonEmpty, SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def output(self):
        return [LocalTarget(os.path.join(get_ext(self.input().path)[0]) + "_1.fastq.gz"),
                LocalTarget(os.path.join(get_ext(self.input().path)[0]) + "_2.fastq.gz")]

    def work_script(self):

        return '''#!/bin/bash
                    source pigz-2.3.3
                    source bedtools-2.17.0;
                    set -euo pipefail

                    bedtools bamtofastq -i {input} \
                                        -fq >(pigz -p 2 -c > {R1}.temp)\
                                        -fq2 >(pigz -p 2 -c > {R2}.temp)

                    mv {R1}.temp {R1}
                    mv {R2}.temp {R2}
        '''.format(input=self.input().path,
                   R1=self.output()[0].path,
                   R2=self.output()[1].path)


@inherits(BAMtoFASTQ, HISAT)
class MergeFQ(CheckTargetNonEmpty, SlurmExecutableTask):
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def requires(self):
        return [self.clone(BAMtoFASTQ, library=lib) for lib in timecourse_libs + PST104E137_rnaseq_libs]

    def input(self):
        logs = [parseHISATLog(self.clone(HISAT, library=lib).output()['log'].path, lib)
                for lib in timecourse_libs + PST104E137_rnaseq_libs]
        filtered = [x['Library'] for x in logs if x['Overall alignment rate'] > 0.1]
        print([(x['Library'], x['Overall alignment rate']) for x in logs if x['Overall alignment rate'] > 10])
        return [self.clone(BAMtoFASTQ, library=lib).output() for lib in filtered]

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, get_ext(os.path.basename(self.reference))[0], 'merged_1.fastq.gz')),
                LocalTarget(os.path.join(self.scratch_dir, get_ext(os.path.basename(self.reference))[0], 'merged_2.fastq.gz'))]

    def work_script(self):
        return '''#!/bin/bash

                    set -euo pipefail

                    cat {input_1} > {R1}.temp
                    cat {input_2} > {R2}.temp

                    mv {R1}.temp {R1}
                    mv {R2}.temp {R2}
        '''.format(input_1=' '.join([x[0].path for x in self.input()]),
                   input_2=' '.join([x[1].path for x in self.input()]),
                   R1=self.output()[0].path,
                   R2=self.output()[1].path)


@requires(MergeFQ)
class Trinity(bioluigi.tasks.trinity.Trinity):
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'trinity.fasta'))


# --------- Portcullis Splicing Junction Filtering --------- #
'''Runs the portcullis pipeline on the merged alignments to produce a
   set of hight quality splice junctions
'''


@requires(MergeBam)
class PortcullisPrep(bioluigi.tasks.portcullis.PortcullisPrep):
    pass


@requires(PortcullisPrep)
class PortcullisJunc(bioluigi.tasks.portcullis.PortcullisJunc):
    pass


@requires(prep=PortcullisPrep, junc=PortcullisJunc)
class PortcullisFilter(bioluigi.tasks.portcullis.PortcullisFilter):
    pass


@requires(ReferenceFasta)
class AddGenome(bioluigi.tasks.jbrowse.AddGenome):
    pass


@requires(genome=AddGenome, portcullis=PortcullisFilter)
class AddPortcullis(bioluigi.tasks.jbrowse.AddPortcullis):
    pass


# --------- Mikado Transcript Merging --------- #
'''Run Mikado pipeline to merge the Stringtie, Cufflinks,
   genome guided Trinity and de novo Trinity transcriptomes
   with the Portcullis splice junctions to make a high confidence
   consensus transcriptome
'''


@requires(cufflinks=AddTranscripts, stringtie=StringTieMerge, portcullis=PortcullisFilter)
class Mikado(bioluigi.tasks.mikado.Mikado):
    pass


@requires(Mikado)
class MikadoGFF3(bioluigi.tasks.misc.GTFtoGFF3):
    pass


@requires(gff=Mikado, genome=AddGenome)
class AddMikadoGFF3(bioluigi.tasks.jbrowse.AddGFF):
    pass


@requires(StringTieGFF3, AddTranscripts, Trinity, TrinityGG, AddPortcullis, AddMikadoGFF3)
class RnaPreprocess(luigi.WrapperTask):
    pass
